/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/InflateError.h"

#include <algorithm>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<InflateError> makerInflateErr_("inflate error");

// -----------------------------------------------------------------------------

void InflateErrorParameters::deserialize(util::CompositePath &path,
                                         const eckit::Configuration &config) {
  oops::Parameters::deserialize(path, config);

  // These checks should really be done at the validation stage (using JSON Schema),
  // but this isn't supported yet, so this is better than nothing.
  if ((inflationFactor.value() == boost::none && inflationVariable.value() == boost::none) ||
      (inflationFactor.value() != boost::none && inflationVariable.value() != boost::none))
    throw eckit::UserError(path.path() +
                           ": Exactly one of the 'inflation factor' and 'inflation variable' "
                           "options must be present");
}

// -----------------------------------------------------------------------------

InflateError::InflateError(const Parameters_ & parameters)
  : allvars_(), parameters_(parameters) {
  if (parameters_.inflationVariable.value() != boost::none) {
    allvars_ += *parameters_.inflationVariable.value();
  }
}

// -----------------------------------------------------------------------------

/// Inflate ObsError by either a constant inflation factor, or by a varying
/// inflation variable.
/// \param vars variables that need to be inflated in the ObsError (filter variables)
/// \param flagged result of "where" statement: which variables/locations need to be
///               updated (has the same variables as \p vars, in the same order)
/// \param data accessor to obs filter data
/// \param flags QC flags (for all "simulated variables")
/// \param obserr ObsError (for all "simulated variables")
void InflateError::apply(const Variables & vars,
                         const std::vector<std::vector<bool>> & flagged,
                         ObsFilterData & data,
                         int /*filterQCflag*/,
                         ioda::ObsDataVector<int> & flags,
                         ioda::ObsDataVector<float> & obserr) const {
  oops::Log::debug() << " InflateError input obserr: " << obserr << std::endl;
  // If float factor is specified
  if (parameters_.inflationFactor.value() != boost::none) {
    float factor = *parameters_.inflationFactor.value();
    for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
      size_t iallvar = obserr.varnames().find(vars.variable(ifiltervar).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[ifiltervar][jobs] && flags[iallvar][jobs] == QCflags::pass) {
          obserr[iallvar][jobs] *= factor;
        }
      }
    }
  // If variable is specified
  } else if (parameters_.inflationVariable.value() != boost::none) {
    const Variable &factorvar = *parameters_.inflationVariable.value();
    ASSERT(factorvar.size() == 1 || factorvar.size() == vars.nvars());
    ioda::ObsDataVector<float> factors(data.obsspace(), factorvar.toOopsObsVariables());
    data.get(factorvar, factors);

    // if inflation factor is 1D variable, apply the same inflation factor to all variables
    // factor_indices = {0, 0, 0, ..., 0} for all nvars
    std::vector<size_t> factor_indices(vars.nvars(), 0);

    // if multiple variables are in the inflation factor, apply different factors to different
    // variables
    // factor_indices = {0, 1, 2, ..., nvars-1}
    if (factorvar.size() == vars.nvars()) {
      std::iota(factor_indices.begin(), factor_indices.end(), 0);
    }

    // loop over all variables to update
    for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
      // find current variable index in obserr
      size_t iallvar = obserr.varnames().find(vars.variable(ifiltervar).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[ifiltervar][jobs] && flags[iallvar][jobs] == QCflags::pass) {
          obserr[iallvar][jobs] *= factors[factor_indices[ifiltervar]][jobs];
        }
      }
    }
  }
  oops::Log::debug() << " InflateError output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
