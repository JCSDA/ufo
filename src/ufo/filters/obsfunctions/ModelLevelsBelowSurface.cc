/*
 * (C) Copyright 2026 NOAA/NCEP/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ModelLevelsBelowSurface.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<ModelLevelsBelowSurface>
  makerModelLevelsBelowSurface_("ModelLevelsBelowSurface");

// -----------------------------------------------------------------------------

ModelLevelsBelowSurface::ModelLevelsBelowSurface(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "ModelLevelsBelowSurface constructor" << std::endl;
  options_.validateAndDeserialize(conf);

  invars_ += options_.obsPressure.value();
  invars_ += Variable("GeoVaLs/air_pressure_at_surface");
  invars_ += Variable("GeoVaLs/air_pressure");
}

// -----------------------------------------------------------------------------

ModelLevelsBelowSurface::~ModelLevelsBelowSurface() {
  oops::Log::trace() << "ModelLevelsBelowSurface destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelLevelsBelowSurface::compute(const ObsFilterData & in,
                                      ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ModelLevelsBelowSurface compute start" << std::endl;

  const float missing = util::missingValue<float>();

  const size_t nlocs = in.nlocs();
  const Variable modelPressure_name("GeoVaLs/air_pressure");
  const size_t nlevs = in.nlevs(modelPressure_name);
  ASSERT_MSG(nlevs >= 2, "ModelLevelsBelowSurface: GeoVaLs/air_pressure must have at least "
                         "2 levels");

  std::vector<float> obsPressure(nlocs);
  in.get(options_.obsPressure.value(), obsPressure);
  std::vector<float> modelSfcPressure(nlocs);
  in.get(Variable("GeoVaLs/air_pressure_at_surface"), modelSfcPressure);

  // get depth of the lowest model layer.
  std::vector<float> modelPrs_lowest(nlocs), modelPrs_second_lowest(nlocs);
  in.get(modelPressure_name, nlevs - 1, modelPrs_lowest);
  in.get(modelPressure_name, nlevs - 2, modelPrs_second_lowest);

  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (obsPressure[jloc] == missing || modelSfcPressure[jloc] == missing ||
        modelPrs_lowest[jloc] == missing || modelPrs_second_lowest[jloc] == missing ||
        modelPrs_lowest[jloc] <= 0.0f || modelPrs_second_lowest[jloc] <= 0.0f) {
      out[0][jloc] = missing;
    } else if (obsPressure[jloc] > modelSfcPressure[jloc]) {
      out[0][jloc] = std::log(obsPressure[jloc] / modelSfcPressure[jloc]) /
                   std::log(modelPrs_lowest[jloc] / modelPrs_second_lowest[jloc]);
    } else {
      out[0][jloc] = 0.0f;
    }
  }
  oops::Log::trace() << "ModelLevelsBelowSurface compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelLevelsBelowSurface::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
