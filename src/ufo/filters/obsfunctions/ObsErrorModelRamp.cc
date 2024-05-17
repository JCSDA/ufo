/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorModelRamp.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

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
  if (options_.x2.value() != boost::none && options_.err2.value() != boost::none) {
    const std::vector<float> &x2 = options_.x2.value().get();
    const std::vector<float> &err2 = options_.err2.value().get();
    ASSERT(x2.size() == yvar.size());
    ASSERT(err2.size() == yvar.size());
    for (size_t i = 0; i < yvar.size(); ++i) {
      ASSERT(x2[i] >= x1[i]);
      ASSERT(err2[i] > 0.0);
    }
  }
}

// -----------------------------------------------------------------------------

ObsErrorModelRamp::~ObsErrorModelRamp() {}

// -----------------------------------------------------------------------------

void ObsErrorModelRamp::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue<float>();

  // Get piece-wise parameters from options
  const std::vector<float> &x0 = options_.x0.value();
  const std::vector<float> &x1 = options_.x1.value();
  const std::vector<float> &err0 = options_.err0.value();
  const std::vector<float> &err1 = options_.err1.value();

  // Check out size
  ASSERT(out.nvars() == x0.size());

  // Compute x values
  const Variable &xvar = options_.xvar.value();
  ioda::ObsDataVector<float> xvals(in.obsspace(), xvar.toOopsObsVariables());
  in.get(xvar, xvals);

  // Optional save of the xfunc values
  if (options_.save) xvals.save("ObsFunction");

  float slope;

  // Optional extra ramp function
  std::vector<float> x2;
  std::vector<float> err2;
  bool cal_err2_x2 = false;
  float slope2 = missing;
  if (options_.x2.value() != boost::none && options_.err2.value() != boost::none) {
    x2 = options_.x2.value().get();
    err2 = options_.err2.value().get();
    cal_err2_x2 = true;
  }
  // Loop over selected variables
  for (size_t jvar = 0; jvar < out.nvars(); ++jvar) {
    size_t ivar = std::min(jvar, xvar.size() - 1);

    // Calculate slope of ramp for this variable
    if (x1[jvar] > x0[jvar]) {
      slope = (err1[jvar] - err0[jvar]) / (x1[jvar] - x0[jvar]);
    } else {
      slope = missing;
    }
    if (cal_err2_x2) {
      if (x2[jvar] > x1[jvar]) {
        slope2 = (err2[jvar] - err1[jvar]) / (x2[jvar] - x1[jvar]);
      } else {
        slope2 = missing;
      }
    }

    // Calculate piece-wise function value across locations
    for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
      out[jvar][iloc] = missing;
      if (xvals[ivar][iloc] != missing) {
        if (xvals[ivar][iloc] <= x0[jvar]) {
          out[jvar][iloc] = err0[jvar];
        } else if (xvals[ivar][iloc] < x1[jvar] && slope != missing) {
          out[jvar][iloc] = err0[jvar] + slope * (xvals[ivar][iloc] - x0[jvar]);
        } else if (cal_err2_x2) {
            if (xvals[ivar][iloc] < x2[jvar] && slope2 != missing) {
              out[jvar][iloc] = err1[jvar] + slope2 * (xvals[ivar][iloc] - x1[jvar]);
            } else {
              out[jvar][iloc] = err2[jvar];
            }
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
