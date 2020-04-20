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
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<AssignError> makerAssignErr_("assign error");

// -----------------------------------------------------------------------------

AssignError::AssignError(const eckit::Configuration & conf)
  : allvars_(), conf_(conf) {
  if (conf_.has("error function")) {
    allvars_ += Variable(conf_.getSubConfiguration("error function"));
  }
  ASSERT(conf_.has("error function") || conf_.has("error parameter"));
}

// -----------------------------------------------------------------------------

void AssignError::apply(const Variables & vars,
                         const std::vector<std::vector<bool>> & flagged,
                         const ObsFilterData & data,
                         ioda::ObsDataVector<int> &,
                         ioda::ObsDataVector<float> & obserr) const {
  oops::Log::debug() << " AssignError input obserr: " << obserr << std::endl;
  // If float error is specified
  if (conf_.has("error parameter")) {
    float error = conf_.getFloat("error parameter");
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[iv][jobs]) obserr[iv][jobs] = error;
      }
    }
  // If variable is specified
  } else if (conf_.has("error function")) {
    Variable errorvar(conf_.getSubConfiguration("error function"));
    ASSERT(errorvar.size() == 1 || errorvar.size() == vars.nvars());
    ioda::ObsDataVector<float> errors(data.obsspace(), errorvar.toOopsVariables(),
                                       errorvar.group(), false);
    data.get(errorvar, errors);
    const float missing = util::missingValue(missing);

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
        if (flagged[iv][jobs] && errors[error_jv[jv]][jobs] != missing)
          obserr[iv][jobs] = errors[error_jv[jv]][jobs];
      }
    }
  }
  oops::Log::debug() << " AssignError output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
