/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ModelLevelIndex.h"

#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

static ObsFunctionMaker<ModelLevelIndex> makerModelLevelIndex_("ModelLevelIndex");

// -----------------------------------------------------------------------------

ModelLevelIndex::ModelLevelIndex(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Add the model vertical coordinate GeoVaL to the list of required variables.
  invars_ += Variable("GeoVaLs/" + options_.modelCoordName.value());
}

// -----------------------------------------------------------------------------

ModelLevelIndex::~ModelLevelIndex() {}

// -----------------------------------------------------------------------------

void ModelLevelIndex::compute(const ObsFilterData & in,
                              ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "ModelLevelIndex::compute started" << std::endl;

  const int missingInt = util::missingValue<int>();
  const float missingFloat = util::missingValue<float>();

  const ioda::ObsSpace & obsdb = in.obsspace();
  const GeoVaLs * const gv = in.getGeoVaLs();

  const std::size_t nlocs = obsdb.nlocs();
  const std::size_t nlevs = gv->nlevs(oops::Variable{options_.modelCoordName.value()});

  // Get observed vertical coordinate.
  std::vector<float> obsVertCoord(nlocs);
  obsdb.get_db(options_.obsCoordGroup, options_.obsCoordName.value(), obsVertCoord);

  std::vector<double> modelVertCoord(nlevs);
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const float z_ob = obsVertCoord[jloc];
    if (z_ob == missingFloat) {
      out[0][jloc] = missingInt;
      continue;
    }

    gv->getAtLocation(modelVertCoord,
                      oops::Variable{options_.modelCoordName.value()},
                      jloc);

    // Out of bounds values are set to missing.
    if (z_ob < *std::min_element(modelVertCoord.begin(), modelVertCoord.end()) ||
        z_ob > *std::max_element(modelVertCoord.begin(), modelVertCoord.end())) {
      out[0][jloc] = missingInt;
      continue;
    }

    const auto[idx, weight] =
      ufo::PiecewiseLinearInterpolation::interpolationIndexAndWeight(modelVertCoord, z_ob);

    out[0][jloc] = idx;
  }

  oops::Log::trace() << "ModelLevelIndex::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelLevelIndex::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
