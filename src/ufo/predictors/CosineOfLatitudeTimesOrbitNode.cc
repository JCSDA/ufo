/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/CosineOfLatitudeTimesOrbitNode.h"

#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<CosineOfLatitudeTimesOrbitNode>
       makerFuncCosineOfLatitudeTimesOrbitNode_("cosine_of_latitude_times_orbit_node");

// -----------------------------------------------------------------------------

CosineOfLatitudeTimesOrbitNode::CosineOfLatitudeTimesOrbitNode(const Parameters_ & parameters,
                                                               const oops::Variables & vars)
  : PredictorBase(parameters, vars) {
}

// -----------------------------------------------------------------------------
void CosineOfLatitudeTimesOrbitNode::compute(const ioda::ObsSpace & odb,
                                             const GeoVaLs &,
                                             const ObsDiagnostics &,
                                             ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();

  // retrieve the sensor view angle
  std::vector<float> cenlat(nlocs, 0.0);
  std::vector<float> node(nlocs, 0.0);
  odb.get_db("MetaData", "latitude", cenlat);
  odb.get_db("MetaData", "sensor_azimuth_angle", node);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = node[jloc] * cos(cenlat[jloc] * Constants::deg2rad);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
