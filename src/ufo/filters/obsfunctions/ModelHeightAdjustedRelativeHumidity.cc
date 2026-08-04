/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ModelHeightAdjustedRelativeHumidity.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ModelHeightAdjustedRelativeHumidity>
        makerModelHeightAdjustedRelativeHumidity_("ModelHeightAdjustedRelativeHumidity");

// -----------------------------------------------------------------------------

ModelHeightAdjustedRelativeHumidity::ModelHeightAdjustedRelativeHumidity(
        const eckit::LocalConfiguration & conf): invars_() {
  oops::Log::trace() << "ModelHeightAdjustedRelativeHumidity constructor" << std::endl;
  // Retrieve observation data  // Required observation data
  invars_ += Variable("ObsValue/relativeHumidityAt2M");
  // Required model surface altitude
  invars_ += Variable("GeoVaLs/height_above_mean_sea_level_at_surface");

  // Required observation station height and temperature
  parameters_.validateAndDeserialize(conf);
  observationRelativeHumidityAtSaturation_ = (
    parameters_.observation_relative_humidity_units.value() ==
    RelativeHumidityUnits::PERCENTAGE) ? 100.0f : 1.0f;
  const Variable elevation = parameters_.elevation.value();
  const Variable temperature = parameters_.temperature.value();
  invars_ += elevation;
  invars_ += temperature;
}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedRelativeHumidity::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ModelHeightAdjustedRelativeHumidity compute start" << std::endl;
  const size_t nlocs = in.nlocs();
  std::vector<float> rh(nlocs);
  std::vector<float> T(nlocs);
  std::vector<float> ModelHeight(nlocs);
  std::vector<float> StationHeight(nlocs);

  in.get(Variable("ObsValue/relativeHumidityAt2M"), rh);
  in.get(parameters_.temperature.value(), T);
  in.get(Variable("GeoVaLs/height_above_mean_sea_level_at_surface"), ModelHeight);
  in.get(parameters_.elevation.value(), StationHeight);

  const float missing = util::missingValue<float>();

  // Maximum values of RH_ice as a fraction for temperatures 0 to -40 deg C
  std::vector<float> rhmax{1.0000, 1.0098, 1.0197, 1.0296, 1.0397, 1.0499,
                           1.0601, 1.0705, 1.0810, 1.0916, 1.1023, 1.1131,
                           1.1240, 1.1351, 1.1462, 1.1575, 1.1688, 1.1803,
                           1.1919, 1.2036, 1.2154, 1.2274, 1.2394, 1.2515,
                           1.2638, 1.2762, 1.2887, 1.3012, 1.3139, 1.3267,
                           1.3396, 1.3526, 1.3658, 1.3790, 1.3923, 1.4057,
                           1.4192, 1.4327, 1.4464, 1.4602, 1.4740};

  // Convert rhmax to percentage if observation variable is a percentage.
  std::transform(rhmax.begin(), rhmax.end(), rhmax.begin(),
            [this](double rh) -> double {return rh * observationRelativeHumidityAtSaturation_;});

  // compute relative humidity correction and adjusted relative humidity.
  const float rhCorrectionPerMeter = 0.0001f * observationRelativeHumidityAtSaturation_;
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (StationHeight[jj] == missing || ModelHeight[jj] == missing || rh[jj] == missing) {
      out[0][jj] = missing;
    } else {
      int Tbin = std::ceil(Constants::t0c - T[jj]);
      Tbin = std::max(0, std::min(40, Tbin));
      float CorrectedRH = std::max(
        rh[jj] - rhCorrectionPerMeter*(StationHeight[jj] - ModelHeight[jj]), 0.0f);
      CorrectedRH = std::min(CorrectedRH, rhmax[Tbin]);
      out[0][jj] = CorrectedRH;
    }
  }
  oops::Log::trace() << "ModelHeightAdjustedRelativeHumidity compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedRelativeHumidity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
