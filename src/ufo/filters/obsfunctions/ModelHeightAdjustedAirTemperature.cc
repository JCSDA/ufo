/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ModelHeightAdjustedAirTemperature.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ModelHeightAdjustedAirTemperature>
        makerModelHeightAdjustedAirTemperature_("ModelHeightAdjustedAirTemperature");

// -----------------------------------------------------------------------------

ModelHeightAdjustedAirTemperature::ModelHeightAdjustedAirTemperature(
        const eckit::LocalConfiguration & conf): invars_() {
  oops::Log::trace() << "ModelHeightAdjustedAirTemperature constructor" << std::endl;

  // Validate and deserialize parameters from YAML config
  parameters_.validateAndDeserialize(conf);

  // Required observation data (configurable, defaults to ObsValue/airTemperatureAt2M)
  std::string obsVarFullName = "ObsValue/airTemperatureAt2M";
  if (parameters_.observedVariable.value() != boost::none) {
    obsVarFullName = parameters_.observedVariable.value()->name.value();
  }
  invars_ += Variable(obsVarFullName);

  // Required model height (terrain or lowest model level)
  std::string modelHeightFrom = "terrain";
  if (parameters_.modelHeight.value() != boost::none) {
    modelHeightFrom = parameters_.modelHeight.value()->from.value();
  }
  if (modelHeightFrom == "terrain") {
    invars_ += Variable("GeoVaLs/height_above_mean_sea_level_at_surface");
  } else if (modelHeightFrom == "lowestModelLevel") {
    invars_ += Variable("GeoVaLs/height_above_mean_sea_level");
  } else {
    throw eckit::BadValue("ModelHeightAdjustedAirTemperature: unknown model height 'from' value '" +
                     modelHeightFrom + "'. Must be 'terrain' or 'lowestModelLevel'.", Here());
  }

  // Required observation station height
  const Variable elevation = parameters_.elevation.value();
  invars_ += elevation;

  // If using local lapse rate, also need model temperature and height profiles
  if (parameters_.lapseRate.value() != boost::none) {
    const std::string scheme = parameters_.lapseRate.value()->scheme.value();
    if (scheme == "local") {
      invars_ += Variable("GeoVaLs/air_temperature");
      // if using the lowestModelLevel, we already have the height above mean sea level;
      // so only add it if otherwise
      if (modelHeightFrom != "lowestModelLevel") {
        invars_ += Variable("GeoVaLs/height_above_mean_sea_level");
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedAirTemperature::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ModelHeightAdjustedAirTemperature compute start" << std::endl;
  const size_t nlocs = in.nlocs();

  // Determine observed variable name
  std::string obsVarFullName = "ObsValue/airTemperatureAt2M";
  if (parameters_.observedVariable.value() != boost::none) {
    obsVarFullName = parameters_.observedVariable.value()->name.value();
  }

  // Determine model height option and offset
  std::string modelHeightFrom = "terrain";
  float modelHeightOffset = 0.0f;
  if (parameters_.modelHeight.value() != boost::none) {
    modelHeightFrom = parameters_.modelHeight.value()->from.value();
    modelHeightOffset = parameters_.modelHeight.value()->offsetInMeters.value();
  }

  // Determine lapse rate scheme
  std::string scheme = "constant";
  if (parameters_.lapseRate.value() != boost::none) {
    scheme = parameters_.lapseRate.value()->scheme.value();
  }

  // Read common input data
  std::vector<float> obsValue(nlocs);
  std::vector<float> ModelHeight(nlocs);
  std::vector<float> StationHeight(nlocs);

  in.get(Variable(obsVarFullName), obsValue);
  if (modelHeightFrom == "terrain") {
    in.get(Variable("GeoVaLs/height_above_mean_sea_level_at_surface"), ModelHeight);
  } else {
    // lowestModelLevel: get height at the lowest model level
    const size_t nlevs = in.nlevs(Variable("GeoVaLs/height_above_mean_sea_level"));
    in.get(Variable("GeoVaLs/height_above_mean_sea_level"), nlevs - 1, ModelHeight);
  }
  in.get(parameters_.elevation.value(), StationHeight);

  // Apply offset to model height, if value is provided
  const float missing = util::missingValue<float>();
  if (modelHeightOffset != 0.0f) {
    for (size_t jj = 0; jj < nlocs; ++jj) {
      if (ModelHeight[jj] != missing) {
        ModelHeight[jj] += modelHeightOffset;
      }
    }
  }

  if (scheme == "constant") {
    // Constant lapse rate mode, use configured value (unit: K/km) if provided,
    // otherwise use Constants::Lclr (unit: K/m)
    float lapseRate;
    if (parameters_.lapseRate.value() != boost::none) {
      lapseRate = parameters_.lapseRate.value()->value.value() / 1000.0f;
    } else {
      lapseRate = static_cast<float>(Constants::Lclr);
    }

    for (size_t jj = 0; jj < nlocs; ++jj) {
      if (StationHeight[jj] == missing || ModelHeight[jj] == missing ||
          obsValue[jj] == missing) {
        out[0][jj] = missing;
      } else {
        out[0][jj] = obsValue[jj] + lapseRate * (StationHeight[jj] - ModelHeight[jj]);
      }
    }
  } else if (scheme == "local") {
    // Local lapse rate mode: compute from model temperature profile
    const int level1 = parameters_.lapseRate.value()->level1.value();
    const int level2 = parameters_.lapseRate.value()->level2.value();

    // Check if bounds are specified for clamping
    const bool useBounds =
        (parameters_.lapseRate.value()->bounds.value() != boost::none);
    float minBound = 0.0f;
    float maxBound = 0.0f;
    if (useBounds) {
      // Convert from K/km to K/m
      minBound = parameters_.lapseRate.value()->bounds.value()->minvalue.value() / 1000.0f;
      maxBound = parameters_.lapseRate.value()->bounds.value()->maxvalue.value() / 1000.0f;
    }

    // Get number of model levels
    const size_t nlevs = in.nlevs(Variable("GeoVaLs/air_temperature"));

    // Level indices counted from surface (1-based) to GeoVaLs indices
    // GeoVaLs convention: level 0 = top of atmosphere, nlevs-1 = surface
    const size_t geoIdx1 = nlevs - static_cast<size_t>(level1);
    const size_t geoIdx2 = nlevs - static_cast<size_t>(level2);

    // Get temperature and height at the two specified levels
    std::vector<float> T1(nlocs), T2(nlocs);
    std::vector<float> Z1(nlocs), Z2(nlocs);

    in.get(Variable("GeoVaLs/air_temperature"), geoIdx1, T1);
    in.get(Variable("GeoVaLs/air_temperature"), geoIdx2, T2);
    in.get(Variable("GeoVaLs/height_above_mean_sea_level"), geoIdx1, Z1);
    in.get(Variable("GeoVaLs/height_above_mean_sea_level"), geoIdx2, Z2);

    for (size_t jj = 0; jj < nlocs; ++jj) {
      if (StationHeight[jj] == missing || ModelHeight[jj] == missing ||
          obsValue[jj] == missing || T1[jj] == missing || T2[jj] == missing ||
          Z1[jj] == missing || Z2[jj] == missing) {
        out[0][jj] = missing;
      } else {
        float dZ = Z2[jj] - Z1[jj];
        float localLapseRate;
        if (std::abs(dZ) < 1.0f) {
          // Safeguard: height difference too small, fall back to standard lapse rate
          localLapseRate = static_cast<float>(Constants::Lclr);
        } else {
          // Lapse rate = (T_lower - T_upper) / (Z_upper - Z_lower) in K/m
          // level1 is closer to surface (lower), level2 is higher
          localLapseRate = (T1[jj] - T2[jj]) / dZ;
        }

        // Apply bounds clamping if specified
        if (useBounds) {
          localLapseRate = std::max(minBound, std::min(maxBound, localLapseRate));
        }

        out[0][jj] = obsValue[jj] + localLapseRate * (StationHeight[jj] - ModelHeight[jj]);
      }
    }
  } else {
    throw eckit::BadValue("ModelHeightAdjustedAirTemperature: unknown lapse rate scheme '" +
                  scheme + "'. Must be 'constant' or 'local'.", Here());
  }

  oops::Log::trace() << "ModelHeightAdjustedAirTemperature compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedAirTemperature::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
