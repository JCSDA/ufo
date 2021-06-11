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

namespace ufo {

static ObsFunctionMaker<LinearCombination>
                       makerLinearCombination_("LinearCombination");

// -----------------------------------------------------------------------------

LinearCombination::LinearCombination(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // number of input variables
  int nv = options_.variables.value().size();

  // get variable informations
  std::vector<ufo::Variable> variables = options_.variables.value();
  for (size_t ii = 0; ii < nv; ++ii) {
    invars_ += variables[ii];
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
  int nv = invars_.size();

  // get coefs for linear combination
  const std::vector<float> coefs = options_.coefs.value();

  // sanity check
  ASSERT(coefs.size() == invars_.size() );

  // compute linear conbination of input variables
  const float missing = util::missingValue(missing);
  std::vector<float> varin(nlocs);
  for (size_t ii = 0; ii < nv; ++ii) {
    in.get(invars_[ii], varin);
    for (size_t jj = 0; jj < nlocs; ++jj) {
      if ( varin[jj] == missing || out[0][jj] == missing ) {
        out[0][jj] = missing;
      } else {
        out[0][jj] += coefs[ii]*varin[jj];
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & LinearCombination::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
