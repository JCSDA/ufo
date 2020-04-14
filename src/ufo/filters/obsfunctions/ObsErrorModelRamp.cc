/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorModelRamp.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorModelRamp>
         makerRamp_("ObsErrorModelRamp");

// -----------------------------------------------------------------------------

ObsErrorModelRamp::ObsErrorModelRamp(const eckit::LocalConfiguration config)
  : invars_() {
  // Initialize options
  options_.deserialize(config);

  // Initialize y-variable
  eckit::LocalConfiguration yconf;
  yconf.set("name", classname());
  if (options_.chlist.value() != boost::none) {
    yconf.set("channels", options_.chlist.value().get());
  }
  yconf.set("options", config);
  Variable yvar(yconf);

  // Initialize x-variable
  const Variable &xvar = options_.xvar.value();
  ASSERT(xvar.size() == 1 ||
         xvar.size() == yvar.size());
  invars_ += xvar;

  // Get piece-wise parameters from options
  const std::vector<float> &x0 = options_.x0.value();
  const std::vector<float> &x1 = options_.x1.value();
  const std::vector<float> &err0 = options_.err0.value();
  const std::vector<float> &err1 = options_.err1.value();

  // Check parameter sizes/values
  ASSERT(yvar.size() == x0.size());
  ASSERT(yvar.size() == x1.size());
  ASSERT(yvar.size() == err0.size());
  ASSERT(yvar.size() == err1.size());
  for (size_t i = 0; i < yvar.size(); ++i) {
    ASSERT(x1[i] >= x0[i]);
    ASSERT(err0[i] > 0.0);
    ASSERT(err1[i] > 0.0);
  }
}

// -----------------------------------------------------------------------------

ObsErrorModelRamp::~ObsErrorModelRamp() {}

// -----------------------------------------------------------------------------

void ObsErrorModelRamp::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue(missing);

  // Get piece-wise parameters from options
  const std::vector<float> &x0 = options_.x0.value();
  const std::vector<float> &x1 = options_.x1.value();
  const std::vector<float> &err0 = options_.err0.value();
  const std::vector<float> &err1 = options_.err1.value();

  // Check out size
  ASSERT(out.nvars() == x0.size());

  // Compute x values
  const Variable &xvar = options_.xvar.value();
  ioda::ObsDataVector<float> xvals(in.obsspace(), xvar.toOopsVariables());
  in.get(xvar, xvals);

  // Optional save of the xfunc values
  if (options_.save) xvals.save("ObsFunction");

  float slope;

  // Loop over selected variables
  for (size_t jvar = 0; jvar < out.nvars(); ++jvar) {
    size_t ivar = std::min(jvar, xvar.size() - 1);

    // Calculate slope of ramp for this variable
    if (x1[jvar] > x0[jvar]) {
      slope = (err1[jvar] - err0[jvar]) / (x1[jvar] - x0[jvar]);
    } else {
      slope = missing;
    }

    // Calculate piece-wise function value across locations
    for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
      out[jvar][iloc] = missing;
      if (xvals[ivar][iloc] != missing) {
        if (xvals[ivar][iloc] <= x0[jvar]) {
          out[jvar][iloc] = err0[jvar];
        } else if (xvals[ivar][iloc] < x1[jvar] && slope != missing) {
          out[jvar][iloc] = err0[jvar] + slope * (xvals[ivar][iloc] - x0[jvar]);
        } else {
          out[jvar][iloc] = err1[jvar];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorModelRamp::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
