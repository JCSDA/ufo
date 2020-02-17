/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/InflateError.h"

#include <algorithm>
#include <set>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<InflateError> makerInflateErr_("inflate error");

// -----------------------------------------------------------------------------

InflateError::InflateError(const eckit::Configuration & conf)
  : allvars_(), conf_(conf) {
  if (conf_.has("inflation variable")) {
    allvars_ += Variable(conf_.getSubConfiguration("inflation variable"));
  }
}

// -----------------------------------------------------------------------------

void InflateError::apply(const Variables & vars,
                         const std::vector<std::vector<bool>> & flagged,
                         const ObsFilterData & data,
                         ioda::ObsDataVector<int> &,
                         ioda::ObsDataVector<float> & obserr) const {
  oops::Log::debug() << " input obserr: " << obserr << std::endl;
  // If float factor is specified
  if (conf_.has("inflation factor")) {
    float factor = conf_.getFloat("inflation factor");
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[iv][jobs]) obserr[iv][jobs] *= factor;
      }
    }
  // If variable is specified
  } else if (conf_.has("inflation variable")) {
    Variable factorvar(conf_.getSubConfiguration("inflation variable"));
    ASSERT(factorvar.size() == 1 || factorvar.size() == vars.nvars());
    oops::Log::debug() << "processing data: " << strfactor_ << std::endl;
    ioda::ObsDataVector<float> factors(data.obsspace(), factorvar.toOopsVariables(),
                                       factorvar.group(), false);
    data.get(factorvar, factors);
    // if inflation factor is 1D variable, apply the same inflation factor to all variables
    // factor_jv = {0, 0, 0, ..., 0} for all nvars
    std::vector<size_t> factor_jv(vars.nvars(), 0);
    // if multiple variables are in the inflation factor, apply different factors to different
    // variables
    // factor_jv = {0, 1, 2, ..., nvars-1}
    if (factorvar.size() == vars.nvars()) {
      std::iota(factor_jv.begin(), factor_jv.end(), 0);
    }
    // loop over all variables to update
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      // find current variable index in obserr
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[iv][jobs]) obserr[iv][jobs] *= factors[factor_jv[jv]][jobs];
      }
    }
  }
  oops::Log::debug() << " output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
