/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/Emissivity_Diff_GMI.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<Emissivity_Diff_GMI> makerEmissivity_Diff_GMI_("Emissivity_Diff_GMI");

Emissivity_Diff_GMI::Emissivity_Diff_GMI(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Check required parameters
  ASSERT(options_.channel.value() != boost::none &&
         options_.regression_constant_1.value() != boost::none &&
         options_.regression_constant_2.value() != boost::none &&
         options_.regression_coeff_1.value() != boost::none &&
         options_.regression_coeff_2.value() != boost::none);

  invars_ += Variable("GeoVaLs/average_surface_temperature_within_field_of_view");
  invars_ += Variable("GeoVaLs/water_area_fraction");
  // GMI channels
  const std::vector<int> channels = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  invars_ += Variable("ObsValue/brightnessTemperature", channels);
}

// -----------------------------------------------------------------------------

Emissivity_Diff_GMI::~Emissivity_Diff_GMI() {}

// -----------------------------------------------------------------------------

void Emissivity_Diff_GMI::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue<float>();
  std::vector<float> &regression_diff = out[0];
  // Get dimension
  const size_t nlocs = in.nlocs();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  // Get average surface temperature in FOV
  std::vector<float> tsavg(nlocs);
  in.get(Variable("GeoVaLs/average_surface_temperature_within_field_of_view"), tsavg);

  const int channel = options_.channel.value().get();
  const float regression_constant_1 = options_.regression_constant_1.value().get();
  const float regression_constant_2 = options_.regression_constant_2.value().get();
  const std::vector<float> &regression_coeff_1 = options_.regression_coeff_1.value().get();
  const std::vector<float> &regression_coeff_2 = options_.regression_coeff_2.value().get();

  std::vector<float> regression_1(nlocs);
  std::vector<float> regression_2(nlocs);
  const std::vector<int> channels = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  std::vector<float> bt_obs(nlocs);
  in.get(Variable("ObsValue/brightnessTemperature", channels), bt_obs);

  // Regression 1
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    regression_1[iloc] = regression_constant_1;
  }
  for (size_t ich = 0; ich < channels.size(); ++ich) {
    in.get(Variable("ObsValue/brightnessTemperature", channels)[ich], bt_obs);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (water_frac[iloc] > 0.99) {
        regression_1[iloc] += bt_obs[iloc] * regression_coeff_1[ich];
      }
    }
  }
  in.get(Variable("ObsValue/brightnessTemperature", channels)[channel-1], bt_obs);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (water_frac[iloc] > 0.99) {
      // Regression 2
      regression_2[iloc] = regression_constant_2;
      regression_2[iloc] += bt_obs[iloc] * regression_coeff_2[0];
      regression_2[iloc] += tsavg[iloc] * regression_coeff_2[1];
      // Difference
      regression_diff[iloc] = regression_1[iloc] - regression_2[iloc];
    } else {
      regression_diff[iloc] = missing;
    }
  }
}
// -----------------------------------------------------------------------------
const ufo::Variables & Emissivity_Diff_GMI::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
