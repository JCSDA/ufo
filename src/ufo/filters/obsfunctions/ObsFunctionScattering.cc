/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionScattering.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionScattering> makerObsFuncScattering_("Scattering");

// -----------------------------------------------------------------------------

ObsFunctionScattering::ObsFunctionScattering() : invars_() {
  // empiracal formula is used to calculate AMSU-A scattering over ocean
  invars_ += "brightness_temperature_1@ObsValue";
  invars_ += "brightness_temperature_2@ObsValue";
  invars_ += "brightness_temperature_15@ObsValue";
}

// -----------------------------------------------------------------------------

ObsFunctionScattering::~ObsFunctionScattering() {}

// -----------------------------------------------------------------------------

void ObsFunctionScattering::compute(const ObsFilterData & input,
                                    ioda::ObsDataVector<float> & out) const {
  // TODO(AS): should use constants for variable names
  const size_t nlocs = input.nlocs();
  std::vector<float> bt1, bt2, bt15;
  input.get("brightness_temperature_1@ObsValue", bt1);
  input.get("brightness_temperature_2@ObsValue", bt2);
  input.get("brightness_temperature_15@ObsValue", bt15);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = -113.2+(2.41-0.0049*bt1[jj])*bt1[jj]+0.454*bt2[jj]-bt15[jj];
    oops::Log::debug() << "Tb1, Tb2, Tb15: " << bt1[jj] << ", " << bt2[jj] << ", " << bt15[jj]
                       << ", scattering=" << out[0][jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionScattering::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
