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
#include "ufo/filters/obsfunctions/ModelHeightAdjustedWindVectorComponent.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// instantiation of the template with northwardWind = false
static ObsFunctionMaker<ModelHeightAdjustedWindVectorComponent<false>>
        eastwardMaker_("ModelHeightAdjustedEastwardWind");
// instantiation of the template with northwardWind = true
static ObsFunctionMaker<ModelHeightAdjustedWindVectorComponent<true>>
        northwardMaker_("ModelHeightAdjustedNorthwardWind");

// -----------------------------------------------------------------------------
template <bool northwardWind>
ModelHeightAdjustedWindVectorComponent<northwardWind>::ModelHeightAdjustedWindVectorComponent(
        const eckit::LocalConfiguration & conf): invars_() {
  // Required observation data
  if (northwardWind) {
    invars_ += Variable("ObsValue/windNorthwardAt10M");
  } else {
    invars_ += Variable("ObsValue/windEastwardAt10M");
  }

  // Required model surface altitude
  invars_ += Variable("GeoVaLs/surface_altitude");

  // Required observation station height
  parameters_.validateAndDeserialize(conf);
  const Variable elevation = parameters_.elevation.value();
  invars_ += elevation;
}

// -----------------------------------------------------------------------------

template <bool northwardWind>
void ModelHeightAdjustedWindVectorComponent<northwardWind>::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  std::vector<float> WindComponent(nlocs);
  std::vector<float> ModelHeight(nlocs);
  std::vector<float> StationHeight(nlocs);

  if (northwardWind) {
    in.get(Variable("ObsValue/windNorthwardAt10M"), WindComponent);
  } else {
    in.get(Variable("ObsValue/windEastwardAt10M"), WindComponent);
  }
  in.get(Variable("GeoVaLs/surface_altitude"), ModelHeight);
  in.get(parameters_.elevation.value(), StationHeight);

  const float missing = util::missingValue<float>();

  // compute wind correction and adjusted winds.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (StationHeight[jj] == missing || ModelHeight[jj] == missing
            || WindComponent[jj] == missing) {
      out[0][jj] = missing;
    } else {
      float HeightDiff;
      HeightDiff = StationHeight[jj] - ModelHeight[jj];
      if (HeightDiff > 100.0) {
        float ScaleFactor;
        ScaleFactor = 1.0/(1.0 + std::min(2.0, (HeightDiff - 100.0) * 0.002));
        out[0][jj] = WindComponent[jj]*ScaleFactor;
      } else {
        out[0][jj] = WindComponent[jj];
      }
    }
  }
}

// -----------------------------------------------------------------------------

template <bool northwardWind>
const ufo::Variables &
    ModelHeightAdjustedWindVectorComponent<northwardWind>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
