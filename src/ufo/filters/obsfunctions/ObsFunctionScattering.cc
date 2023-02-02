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
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionScattering> makerObsFuncScattering_("Scattering");

// -----------------------------------------------------------------------------

ObsFunctionScattering::ObsFunctionScattering(const eckit::LocalConfiguration)
  : invars_() {
  // empiracal formula is used to calculate AMSU-A scattering over ocean
  std::vector<int> channels{1, 2, 15};
  invars_ += Variable("ObsValue/brightnessTemperature", channels);
}

// -----------------------------------------------------------------------------

ObsFunctionScattering::~ObsFunctionScattering() {}

// -----------------------------------------------------------------------------

void ObsFunctionScattering::compute(const ObsFilterData & input,
                                    ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = input.nlocs();
  std::vector<float> bt1, bt2, bt15;
  std::vector<int> channels{1, 2, 15};

  input.get(Variable("ObsValue/brightnessTemperature", channels)[0], bt1);
  input.get(Variable("ObsValue/brightnessTemperature", channels)[1], bt2);
  input.get(Variable("ObsValue/brightnessTemperature", channels)[2], bt15);

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
