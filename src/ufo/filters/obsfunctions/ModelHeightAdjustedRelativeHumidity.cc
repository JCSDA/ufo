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
  // Retrieve observation data  // Required observation data
  invars_ += Variable("ObsValue/relativeHumidityAt2M");
  // Required model surface altitude
  invars_ += Variable("GeoVaLs/surface_altitude");

  // Required observation station height and temperature
  parameters_.validateAndDeserialize(conf);
  const Variable elevation = parameters_.elevation.value();
  const Variable temperature = parameters_.temperature.value();
  invars_ += elevation;
  invars_ += temperature;
}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedRelativeHumidity::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  std::vector<float> rh(nlocs);
  std::vector<float> T(nlocs);
  std::vector<float> ModelHeight(nlocs);
  std::vector<float> StationHeight(nlocs);

  in.get(Variable("ObsValue/relativeHumidityAt2M"), rh);
  in.get(parameters_.temperature.value(), T);
  in.get(Variable("GeoVaLs/surface_altitude"), ModelHeight);
  in.get(parameters_.elevation.value(), StationHeight);

  const float missing = util::missingValue<float>();

  // Maximum values of RH_ice for temperatures 0 to -40 deg C
  const std::vector<float>  rhmax{100.00, 100.98, 101.97, 102.96, 103.97, 104.99,
                                  106.01, 107.05, 108.10, 109.16, 110.23, 111.31,
                                  112.40, 113.51, 114.62, 115.75, 116.88, 118.03,
                                  119.19, 120.36, 121.54, 122.74, 123.94, 125.15,
                                  126.38, 127.62, 128.87, 130.12, 131.39, 132.67,
                                  133.96, 135.26, 136.58, 137.90, 139.23, 140.57,
                                  141.92, 143.27, 144.64, 146.02, 147.40};

  // compute relative humidity correction and adjusted relative humidity.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (StationHeight[jj] == missing || ModelHeight[jj] == missing || rh[jj] == missing) {
      out[0][jj] = missing;
    } else {
      int Tbin = std::ceil(Constants::t0c - T[jj]);
      Tbin = std::max(0, std::min(40, Tbin));
      float CorrectedRH = std::max(rh[jj] - 0.01*(StationHeight[jj] - ModelHeight[jj]), 0.0);
      CorrectedRH = std::min(CorrectedRH, rhmax[Tbin]);
      out[0][jj] = CorrectedRH;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedRelativeHumidity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
