/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionVelocity.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionVelocity> makerObsFuncVelocity_("Velocity");

// -----------------------------------------------------------------------------

ObsFunctionVelocity::ObsFunctionVelocity(const eckit::LocalConfiguration conf)
  : invars_(), group_() {
  group_ = conf.getString("type", "ObsValue");
  invars_ += Variable("eastward_wind@" + group_);
  invars_ += Variable("northward_wind@" + group_);
}

// -----------------------------------------------------------------------------

ObsFunctionVelocity::~ObsFunctionVelocity() {}

// -----------------------------------------------------------------------------

void ObsFunctionVelocity::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue(missing);
  std::vector<float> u, v;
  in.get(Variable("eastward_wind@" + group_), u);
  in.get(Variable("northward_wind@" + group_), v);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      out[0][jj] = sqrt(u[jj]*u[jj] + v[jj]*v[jj]);
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionVelocity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
