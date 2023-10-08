/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Function to calculate the average model temperature below a given height.
 * Can be calculated as a simple average over model levels, or take into account
 * the thickness of each model layer.
 * -----------------------------------------------------------------------------
 */
#include "ufo/filters/obsfunctions/AverageTemperatureBelow.h"

#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<AverageTemperatureBelow>
  makerAverageTemperatureBelow_("AverageTemperatureBelow");

/* -----------------------------------------------------------------------------
 * Specify that air_temperature and geopotential_height need to be
 * provided to this function
 * -----------------------------------------------------------------------------
 */
AverageTemperatureBelow::AverageTemperatureBelow(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Get options from argument
  options_.deserialize(conf);

  invars_ += Variable("GeoVaLs/air_temperature");
  invars_ += Variable("GeoVaLs/geopotential_height");
}

// -----------------------------------------------------------------------------

AverageTemperatureBelow::~AverageTemperatureBelow() {}

/* -----------------------------------------------------------------------------
 * Perform the computation.  Read in the geovals for air_temperature and
 * geopotential_height, and use these to find the average temperature below
 * the given level.
 * -----------------------------------------------------------------------------
 */
void AverageTemperatureBelow::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/air_temperature"));

  oops::Log::debug() << "AvereageTemperatureBelow: nlocs=" << nlocs << " nlevs="
                     << nlevs << std::endl;

  std::vector<float> temperature_sum(nlocs, 0.0f);
  std::vector<float> weights_sum(nlocs, 0);

  // Choose code path based on whether we are using the layer thickness to
  // weight the temperature calculation
  if (options_.useThickness.value()) {
    std::vector<float> prev_temp;
    std::vector<float> prev_height;
    in.get(Variable("GeoVaLs/air_temperature"), 0, prev_temp);
    in.get(Variable("GeoVaLs/geopotential_height"), 0, prev_height);
    for (size_t ilev=1; ilev < nlevs; ++ilev) {
      std::vector<float> air_temp;
      std::vector<float> height;
      in.get(Variable("GeoVaLs/air_temperature"), ilev, air_temp);
      in.get(Variable("GeoVaLs/geopotential_height"), ilev, height);
      for (size_t iloc=0; iloc < nlocs; ++iloc) {
        // Work out if we have part of this interval above the given threshold.
        // If so, make the calculation based on the part of the layer which is
        // below the threshold.
        if (height[iloc] < options_.heightLimit.value()) {
          float heightDiff;
          float boundaryTemp;
          if (prev_height[iloc] < options_.heightLimit.value()) {
            heightDiff = prev_height[iloc] - height[iloc];
            boundaryTemp = prev_temp[iloc];
          } else {
            heightDiff = options_.heightLimit.value() - height[iloc];
            boundaryTemp = air_temp[iloc] + (prev_temp[iloc] - air_temp[iloc]) *
                           heightDiff / (prev_height[iloc] - height[iloc]);
          }
          temperature_sum[iloc] += heightDiff * (air_temp[iloc] + boundaryTemp) / 2.0f;
          weights_sum[iloc] += heightDiff;
        }
      }
      prev_temp = air_temp;
      prev_height = height;
    }
  } else {
    // If not using layer thickness, then just calculate the average over the
    // matching model levels, considering each level in turn.
    for (size_t ilev=0; ilev < nlevs; ++ilev) {
      std::vector<float> air_temp;
      std::vector<float> height;
      in.get(Variable("GeoVaLs/air_temperature"), ilev, air_temp);
      in.get(Variable("GeoVaLs/geopotential_height"), ilev, height);
      for (size_t iloc=0; iloc < nlocs; ++iloc) {
        if (height[iloc] < options_.heightLimit.value()) {
          temperature_sum[iloc] += air_temp[iloc];
          weights_sum[iloc]++;
        }
      }
    }
  }

  // Normalise the sum by the total number of levels, or the total layer depth.
  const float missingFloat = util::missingValue<float>();
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    if (weights_sum[iloc] > 0) {
      out[0][iloc] = temperature_sum[iloc] / weights_sum[iloc];
    } else {
      out[0][iloc] = missingFloat;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & AverageTemperatureBelow::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
