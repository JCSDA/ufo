/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorModelQuad.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/PropertiesOfNVectors.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorModelQuad>
         makerQuad_("ObsErrorModelQuad");

// -----------------------------------------------------------------------------

ObsErrorModelQuad::ObsErrorModelQuad(const eckit::LocalConfiguration & config)
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
  const std::vector<float> &a = options_.a.value();
  const std::vector<float> &b = options_.b.value();
  const std::vector<float> &err0 = options_.err0.value();
  const std::vector<float> &err1 = options_.err1.value();

  // Check parameter sizes/values
  if (!oops::allVectorsSameNonZeroSize(a, b, err0, err1)) {
    oops::Log::warning() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(a, b, err0, err1)
                         << std::endl;
    return;
  }

  for (size_t i = 0; i < yvar.size(); ++i) {
    ASSERT(abs(a[i]) > 0.0);
    ASSERT(err0[i] > 0.0);
    ASSERT(err1[i] > 0.0);
    ASSERT(err1[i] >= err0[i]);
  }
}

// -----------------------------------------------------------------------------

ObsErrorModelQuad::~ObsErrorModelQuad() {}

// -----------------------------------------------------------------------------

void ObsErrorModelQuad::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue<float>();

  // Get piece-wise parameters from options
  const std::vector<float> &a = options_.a.value();
  const std::vector<float> &b = options_.b.value();
  const std::vector<float> &err0 = options_.err0.value();
  const std::vector<float> &err1 = options_.err1.value();

  // Check out size
  ASSERT(out.nvars() == a.size());

  // Compute x values
  const Variable &xvar = options_.xvar.value();
  ioda::ObsDataVector<float> xvals(in.obsspace(), xvar.toOopsObsVariables());
  in.get(xvar, xvals);

  // Optional save of the xfunc values
  if (options_.save) xvals.save("ObsFunction");

  float x0, x1, c;

  for (size_t jvar = 0; jvar < out.nvars(); ++jvar) {
    size_t ivar = std::min(jvar, xvar.size() - 1);

    // Calculate the inflection point x-values (x0, x1)
    // and quadratic apex y-value (c)
    if (a[jvar] < 0.0f) {
      c = err1[jvar];
      x1 = b[jvar];
      x0 = x1 - sqrt((err0[jvar] - c) / a[jvar]);
    } else {
      c = err0[jvar];
      x0 = b[jvar];
      x1 = x0 + sqrt((err1[jvar] - c) / a[jvar]);
    }

    // Calculate piece-wise function value across locations
    for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
      out[jvar][iloc] = missing;
      if (xvals[ivar][iloc] != missing) {
        if (xvals[ivar][iloc] <= x0) {
          out[jvar][iloc] = err0[jvar];
        } else if (xvals[ivar][iloc] < x1) {
          out[jvar][iloc] = a[jvar] * pow(xvals[ivar][iloc] - b[jvar], 2) + c;
        } else {
          out[jvar][iloc] = err1[jvar];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorModelQuad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
