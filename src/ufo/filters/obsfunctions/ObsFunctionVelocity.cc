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

namespace ufo {

static ObsFunctionMaker<ObsFunctionVelocity> makerObsFuncVelocity_("Velocity");

// -----------------------------------------------------------------------------

ObsFunctionVelocity::ObsFunctionVelocity()
  : invars_() {
  invars_ += "eastward_wind@ObsValue";
  invars_ += "northward_wind@ObsValue";
}

// -----------------------------------------------------------------------------

ObsFunctionVelocity::~ObsFunctionVelocity() {}

// -----------------------------------------------------------------------------

void ObsFunctionVelocity::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // TODO(AS): should use constants for variable names
  const size_t nlocs = in.nlocs();
  std::vector<float> u = in.get("eastward_wind@ObsValue");
  std::vector<float> v = in.get("northward_wind@ObsValue");
  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = sqrt(pow(u[jj], 2) + pow(v[jj], 2));
    oops::Log::debug() << "u, v: " << u[jj] << ", "
                       << v[jj] << ", speed=" << out[0][jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionVelocity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
