/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ModelLevelIndex.h"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

static ObsFunctionMaker<ModelLevelIndex> makerModelLevelIndex_("ModelLevelIndex");

// -----------------------------------------------------------------------------

ModelLevelIndex::ModelLevelIndex(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "ModelLevelIndex constructor" << std::endl;
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Add the model vertical coordinate GeoVaL to the list of required variables.
  invars_ += Variable("GeoVaLs/" + options_.modelCoordName.value());
}

// -----------------------------------------------------------------------------

ModelLevelIndex::~ModelLevelIndex() {
  oops::Log::trace() << "ModelLevelIndex destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelLevelIndex::compute(const ObsFilterData & in,
                              ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "ModelLevelIndex compute start" << std::endl;

  const int missingInt = util::missingValue<int>();
  const float missingFloat = util::missingValue<float>();

  const ioda::ObsSpace & obsdb = in.obsspace();
  const GeoVaLs * const gv = in.getGeoVaLs();

  const std::size_t nlocs = obsdb.nlocs();
  const std::size_t nlevs = gv->nlevs(oops::Variable{options_.modelCoordName.value()});
  const std::size_t indexOffset = options_.indexModelLevelsFromOne.value() ? 1 : 0;

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

    auto[idx, weight] =
      ufo::PiecewiseLinearInterpolation::interpolationIndexAndWeight(modelVertCoord, z_ob);
    if (options_.closestModelIndex.value() && weight <= 0.5) {
      // If the observation is closer to the adjacent model level (idx + 1), return
      // that index instead.
      idx += 1;
    }

    // Output the user requested index inverting from top-down to bottom-up if requested.
    out[0][jloc] = missingInt;  // Default to missing in case idx is out of bounds.
    const int nlevsInt = static_cast<int>(nlevs);
    const int indexOffsetInt = static_cast<int>(indexOffset);
    if (idx >= 0 && idx < nlevsInt) {
      if (options_.invertModelIndex.value()) {
        out[0][jloc] = nlevsInt - 1 - idx + indexOffsetInt;
      } else {
        out[0][jloc] = idx + indexOffsetInt;
      }
    }
  }
  oops::Log::trace() << "ModelLevelIndex compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelLevelIndex::requiredVariables() const {
  return invars_;
}

}  // namespace ufo
