/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorModelStepwiseLinear.h"

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PropertiesOfNVectors.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorModelStepwiseLinear> makerSteps_("ObsErrorModelStepwiseLinear");

// -----------------------------------------------------------------------------

ObsErrorModelStepwiseLinear::ObsErrorModelStepwiseLinear(const eckit::LocalConfiguration config)
  : invars_() {
  // Initialize options
  options_.deserialize(config);

  // Get piece-wise parameters from options.
  const std::vector<float> &xvals = options_.xvals.value();
  const std::vector<float> &errors = options_.errors.value();

  // Ensure same size vectors (xvals and errors); Also ensure more than one value in each.
  if (!oops::allVectorsSameNonZeroSize(xvals, errors) ||
                ((xvals.size() <= 1 || errors.size() <= 1))) {
      std::ostringstream errString;
      errString << "At least one vector is the wrong size - either unequal or less than 2."
                << std::endl << "Vector sizes of xvals, errors: "
                << oops::listOfVectorSizes(xvals, errors) << std::endl;
      oops::Log::error() << errString.str();
      throw eckit::BadValue(errString.str());
  }

  // Initialize x-variable
  const Variable &xvar = options_.xvar.value();
  ASSERT(xvar.size() == 1);
  invars_ += xvar;

  // Check that all errors >= 0
  for (size_t i = 0; i < errors.size(); ++i) {
    ASSERT(errors[i] >= 0.0);
  }

  // In order to check for beyond range of xvals, determine if ascending or descending.
  // Also ensure that the entire vector is consistent throughout and no consecutive elements
  // of xval are equal.
  if (xvals.back() < xvals.front()) {
    isAscending_ = false;
  }
  for (size_t kv = 1; kv < xvals.size(); ++kv) {
    if ((xvals[kv] >= xvals[kv-1] && !isAscending_) || (xvals[kv] <= xvals[kv-1] && isAscending_)) {
      std::ostringstream errString;
      errString << "The xvals vector of elements is NOT internally consistent."
                << " It must be entirely either ascending or descending order,"
                << " or two consecutive values are the same." << std::endl;
      oops::Log::error() << errString.str();
      throw eckit::BadValue(errString.str());
    }
  }
  oops::Log::debug() << "ObsErrorModelStepwiseLinear: config (constructor) = "
                     << config << std::endl;
}

// -----------------------------------------------------------------------------

ObsErrorModelStepwiseLinear::~ObsErrorModelStepwiseLinear() {}

// -----------------------------------------------------------------------------

void ObsErrorModelStepwiseLinear::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue(missing);
  // Linearly interpolate from y0 to y1 at xstar between x0 and x1 to arrive at error
  float x0, x1, y0, y1;
  float xstar, error;

  // Get the x-variable name and piece-wise parameters from options
  const Variable &xvar = options_.xvar.value();
  const std::vector<float> &xvals = options_.xvals.value();
  const std::vector<float> &errors = options_.errors.value();
  oops::Log::debug() << "  ObsErrorModelStepwiseLinear, x-variable name: " << xvar.variable()
                     << "  and group: " << xvar.group() << std::endl;

  // Populate the testdata array.  xstar is just the 0..nloc-1 value of testvar[iv]
  // At each nloc, find matching (x0,x1) and (y0,y1) pair for linear interp.
  ioda::ObsDataVector<float> testdata(data.obsspace(), xvar.toOopsVariables());
  data.get(xvar, testdata);

  // The 1st index of data should have size 1 and 2nd index should be size nlocs.
  int iv = 0;
  if (testdata[iv].size() != obserr[iv].size()) {
    std::ostringstream errString;
    errString << "Something is wrong, testdata size not equal obserr size."
              << " Sizes: " << testdata[iv].size() << " and " << obserr[iv].size() << std::endl;
    oops::Log::error() << errString.str();
    throw eckit::BadValue(errString.str());
  }

  for (size_t jobs = 0; jobs < testdata[iv].size(); ++jobs) {
    error = 0.0;
    obserr[iv][jobs] = missing;
    if (testdata[iv][jobs] == missing) {
      continue;
    }
    xstar = testdata[iv][jobs];
    if ((xstar <= xvals[0] && isAscending_) || (xstar >= xvals[0] && !isAscending_)) {
      error = errors[0];
    } else if ((xstar >= xvals.back() && isAscending_)
              || (xstar <= xvals.back() && !isAscending_)) {
      error = errors[errors.size()-1];
    } else {
      for (size_t kv = 1; kv < xvals.size(); ++kv) {
        if ((isAscending_ && (xstar > xvals[kv-1]) && (xstar <= xvals[kv])) ||
              (!isAscending_ && (xstar < xvals[kv-1]) && (xstar >= xvals[kv]))) {
          x0 = xvals[kv-1];
          x1 = xvals[kv];
          y0 = errors[kv-1];
          y1 = errors[kv];
          error = y0 + (xstar-x0)*((y1-y0)/(x1-x0));
          break;
        }
      }
    }
    // TODO(gthompsn):  probably need this next line for when filtervariable is flagged missing
    // if (!flagged_[jv][jobs]) obserr[jv][jobs] = error;
    obserr[iv][jobs] = error;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorModelStepwiseLinear::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
