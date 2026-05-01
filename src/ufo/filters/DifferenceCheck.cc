/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/DifferenceCheck.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

namespace {

// -----------------------------------------------------------------------------

/// Set each element of \p flagged to true if the corresponding element of \p apply is true
/// and the corresponding element of \p testValues or \p refValues is missing value.
/// If the corresponding element of \p testValues and \p refValues is non-missing values
/// then set \p flagged to true if their difference is lying outside [\p minValue, \p maxValue].
void flagWhereDiffOut(const std::vector<float> &refValues, std::vector<float> &testValues,
                       const size_t nlocs, const size_t nvars,
                       const float vmin, const float vmax,
                       const bool minExclusive, const bool maxExclusive,
                       const std::vector<bool> & apply,
                       std::vector<std::vector<bool>> & flagged) {
  const float missing = util::missingValue<float>();
  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // check to see if one of the reference or value is missing
      if (testValues[jobs] == missing || refValues[jobs] == missing) {
        for (size_t jv = 0; jv < nvars; ++jv) flagged[jv][jobs] = true;
      } else {
        // Check if difference is within min/max value range and set flag
        float diff = testValues[jobs] - refValues[jobs];
        for (size_t jv = 0; jv < nvars; ++jv) {
          if (vmin != missing && diff <= vmin) {
            if (diff < vmin || minExclusive) flagged[jv][jobs] = true;
          }
          if (!flagged[jv][jobs] && vmax != missing && diff >= vmax) {
            if (diff > vmax || maxExclusive) flagged[jv][jobs] = true;
          }
        }  // end for loop over vars
      }  // end else for missing value check
    }  // end if apply[jobs]
  }  // end for loop over obs
}
}  // end of namespace

// -----------------------------------------------------------------------------

DifferenceCheck::DifferenceCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 ioda::ObsDataVector<int> & flags,
                                 ioda::ObsDataVector<float> & obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::trace() << "DifferenceCheck constructor" << std::endl;
  allvars_ += parameters_.ref;
  allvars_ += parameters_.val;
}

// -----------------------------------------------------------------------------

DifferenceCheck::~DifferenceCheck() {
  oops::Log::trace() << "DifferenceCheck destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void DifferenceCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "DifferenceCheck applyFilter start" << std::endl;

  const float missing = util::missingValue<float>();
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = filtervars.nvars();

// min/max value setup
  float vmin = parameters_.minvalue.value().value_or(missing);
  float vmax = parameters_.maxvalue.value().value_or(missing);
  const bool minExclusive = parameters_.minExclusive.value();
  const bool maxExclusive = parameters_.maxExclusive.value();

// check for threshold and if exists, set vmin and vmax appropriately
  if (parameters_.threshold.value() != boost::none) {
    const float thresh = *parameters_.threshold.value();
    vmin = -thresh;
    vmax = thresh;
  }

  // Get reference values and values as ObsDataVector using variable name with channel
  // to compare (as floats)
  ioda::ObsDataVector<float> refDataVec(obsdb_, parameters_.ref.value().toOopsObsVariables(),
                                         parameters_.ref.value().group(), false);
  ioda::ObsDataVector<float> valDataVec(obsdb_, parameters_.val.value().toOopsObsVariables(),
                                         parameters_.val.value().group(), false);
  data_.get(parameters_.ref, refDataVec);
  data_.get(parameters_.val, valDataVec);
  ASSERT(refDataVec.nlocs() == valDataVec.nlocs());
  ASSERT(refDataVec.nvars() == valDataVec.nvars());

  for (size_t ivar = 0; ivar < refDataVec.nvars(); ++ivar) {
    flagWhereDiffOut(refDataVec[ivar], valDataVec[ivar],
             nlocs, nvars, vmin, vmax, minExclusive, maxExclusive, apply, flagged);
  }

  oops::Log::trace() << "DifferenceCheck applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void DifferenceCheck::print(std::ostream & os) const {
  os << "DifferenceCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
