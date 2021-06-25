/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ModelHeightAdjustedAirTemperature.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ModelHeightAdjustedAirTemperature>
        makerModelHeightAdjustedAirTemperature_("ModelHeightAdjustedAirTemperature");

// -----------------------------------------------------------------------------

ModelHeightAdjustedAirTemperature::ModelHeightAdjustedAirTemperature(
        const eckit::LocalConfiguration & conf): invars_() {
  // Required observation data
  invars_ += Variable("air_temperature_at_2m@ObsValue");
  // Required model surface altitude
  invars_ += Variable("surface_altitude@GeoVaLs");

  // Required observation station height
  parameters_.validateAndDeserialize(conf);
  const Variable elevation = parameters_.elevation.value();
  invars_ += elevation;
}

// -----------------------------------------------------------------------------

ModelHeightAdjustedAirTemperature::~ModelHeightAdjustedAirTemperature() {}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedAirTemperature::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  std::vector<float> t2(nlocs);
  std::vector<float> ModelHeight(nlocs);
  std::vector<float> StationHeight(nlocs);

  in.get(Variable("air_temperature_at_2m@ObsValue"), t2);
  in.get(Variable("surface_altitude@GeoVaLs"), ModelHeight);
  in.get(parameters_.elevation.value(), StationHeight);

  const float missing = util::missingValue(missing);

  // compute temperature correction and adjusted temperature.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (StationHeight[jj] == missing || ModelHeight[jj] == missing || t2[jj] == missing) {
      out[0][jj] = missing;
    } else {
      out[0][jj] = t2[jj] + Constants::Lclr*(StationHeight[jj] - ModelHeight[jj]);
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedAirTemperature::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
