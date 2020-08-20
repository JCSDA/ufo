/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/CosineOfLatitudeTimesOrbitNode.h"

#include "ioda/ObsSpace.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<CosineOfLatitudeTimesOrbitNode>
       makerFuncCosineOfLatitudeTimesOrbitNode_("cosine_of_latitude_times_orbit_node");

// -----------------------------------------------------------------------------

CosineOfLatitudeTimesOrbitNode::CosineOfLatitudeTimesOrbitNode(
                                const eckit::Configuration & conf, const std::vector<int> & jobs)
  : PredictorBase(conf, jobs) {
}

// -----------------------------------------------------------------------------
void CosineOfLatitudeTimesOrbitNode::compute(const ioda::ObsSpace & odb,
                                             const GeoVaLs &,
                                             const ObsDiagnostics &,
                                             ioda::ObsDataVector<double> & out) const {
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  // Following variables should be moved to yaml file ?
  const double ssmis_precond = 0.01;  //  default preconditioner for ssmis bias terms

  // retrieve the sensor view angle
  std::vector<float> cenlat(nlocs, 0.0);
  std::vector<float> node(nlocs, 0.0);
  odb.get_db("MetaData", "latitude", cenlat);
  odb.get_db("MetaData", "sensor_azimuth_angle", node);

  for (const auto & job : jobs_) {
    const std::string varname = name() + "_" + std::to_string(job);
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      out[varname][jl] = node[jl]*cos(cenlat[jl]*Constants::deg2rad);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
