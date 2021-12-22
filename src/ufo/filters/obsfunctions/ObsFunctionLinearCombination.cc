/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionLinearCombination.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<LinearCombination<float>> floatMaker("LinearCombination");
static ObsFunctionMaker<LinearCombination<int>> intMaker("LinearCombination");

// -----------------------------------------------------------------------------

template <typename FunctionValue>
LinearCombination<FunctionValue>::LinearCombination(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variables.value()) {
    invars_ += var;
  }
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void LinearCombination<FunctionValue>::compute(const ObsFilterData & in,
                                               ioda::ObsDataVector<FunctionValue> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // get coefs for linear combination
  const std::vector<FunctionValue> coefs = options_.coefs.value();

  // sanity checks / initialize
  ASSERT(coefs.size() == nv);
  out.zero();

  // compute linear combination of input variables
  const FunctionValue missing = util::missingValue(missing);
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<FunctionValue> varin(in.obsspace(), invars_[ivar].toOopsVariables());
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

template <typename FunctionValue>
const ufo::Variables & LinearCombination<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
