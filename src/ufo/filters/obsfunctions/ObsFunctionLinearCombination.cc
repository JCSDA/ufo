/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionLinearCombination.h"

#include <algorithm>
#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"

namespace ufo {

static ObsFunctionMaker<LinearCombination>
                       makerLinearCombination_("LinearCombination");

// -----------------------------------------------------------------------------

LinearCombination::LinearCombination(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variables.value()) {
    invars_ += var;
  }
}

// -----------------------------------------------------------------------------

LinearCombination::~LinearCombination() {}

// -----------------------------------------------------------------------------

void LinearCombination::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // get coefs for linear combination
  const std::vector<float> coefs = options_.coefs.value();

  // sanity checks / initialize
  ASSERT(coefs.size() == nv);
  out.zero();

  // compute linear combination of input variables
  const float missing = util::missingValue(missing);
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<float> varin(in.obsspace(), invars_[ivar].toOopsVariables());
    in.get(invars_[ivar], varin);
    ASSERT(varin.nvars() == out.nvars());
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < out.nvars(); ++ichan) {
        if ( varin[ichan][iloc] == missing || out[ichan][iloc] == missing ) {
          out[ichan][iloc] = missing;
        } else {
          out[ichan][iloc] += coefs[ivar] * varin[ichan][iloc];
        }
      }  // ichan
    }  // nlocs
  }  // nvars
}

// -----------------------------------------------------------------------------

const ufo::Variables & LinearCombination::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
