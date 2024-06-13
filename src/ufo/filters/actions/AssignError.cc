/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/AssignError.h"

#include <algorithm>
#include <set>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<AssignError> makerAssignErr_("assign error");

// -----------------------------------------------------------------------------

void AssignErrorParameters::deserialize(util::CompositePath &path,
                                        const eckit::Configuration &config) {
  oops::Parameters::deserialize(path, config);

  // These checks should really be done at the validation stage (using JSON Schema),
  // but this isn't supported yet, so this is better than nothing.
  if ((errorParameter.value() == boost::none && errorFunction.value() == boost::none
       && errorParameterVector.value() == boost::none) ||
      (errorParameter.value() != boost::none && errorFunction.value() != boost::none
       && errorParameterVector.value() != boost::none))
    throw eckit::UserError(path.path() +
                           ": Exactly one of the 'error parameter' and 'error function' "
                           "options must be present");
}

// -----------------------------------------------------------------------------

AssignError::AssignError(const Parameters_ & parameters)
  : allvars_(), parameters_(parameters) {
  if (parameters_.errorFunction.value() != boost::none) {
    allvars_ += *parameters_.errorFunction.value();
  }
}

// -----------------------------------------------------------------------------

void AssignError::apply(const Variables & vars,
                        const std::vector<std::vector<bool>> &mask,
                        ObsFilterData & data,
                        int /*filterQCflag*/,
                        ioda::ObsDataVector<int> & qcFlags,
                        ioda::ObsDataVector<float> & obserr) const {
  oops::Log::debug() << " AssignError input obserr: " << obserr << std::endl;
  const float missing = util::missingValue<float>();
  // If float error is specified
  if (parameters_.errorParameter.value() != boost::none) {
    float error = *parameters_.errorParameter.value();
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (mask[jv][jobs] && error != missing) {
            obserr[iv][jobs] = error;
        }
      }
    }
    // If variable is specified
  } else if (parameters_.errorParameterVector.value() != boost::none) {
    std::vector<float> errorvector = *parameters_.errorParameterVector.value();
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (mask[jv][jobs] && errorvector[jv] != missing) {
          obserr[iv][jobs] = errorvector[jv];
        }
      }
    }
    // If variable is specified
  } else if (parameters_.errorFunction.value() != boost::none) {
    const Variable &errorvar = *parameters_.errorFunction.value();
    ASSERT(errorvar.size() == 1 || errorvar.size() == vars.nvars());
    ioda::ObsDataVector<float> errors(data.obsspace(), errorvar.toOopsObsVariables(),
                                      errorvar.group(), false);
    data.get(errorvar, errors);

    // if assigned error function is 1D variable, apply the same error to all variables
    // error_jv = {0, 0, 0, ..., 0} for all nvars
    std::vector<size_t> error_jv(vars.nvars(), 0);
    // if multiple variables are in the assigned error function, apply different error to different
    // variables
    // error_jv = {0, 1, 2, ..., nvars-1}
    if (errorvar.size() == vars.nvars()) {
      std::iota(error_jv.begin(), error_jv.end(), 0);
    }
    // loop over all variables to update
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      // find current variable index in obserr
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (mask[jv][jobs] && errors[error_jv[jv]][jobs] != missing) {
           obserr[iv][jobs] = errors[error_jv[jv]][jobs];
        }
      }
    }
  }
  oops::Log::debug() << " AssignError output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
