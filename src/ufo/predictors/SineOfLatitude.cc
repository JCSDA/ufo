/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/SineOfLatitude.h"

#include "ioda/ObsSpace.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<SineOfLatitude>
       makerFuncSineOfLatitude_("sine_of_latitude");

// -----------------------------------------------------------------------------

SineOfLatitude::SineOfLatitude(const eckit::Configuration & conf, const std::vector<int> & jobs)
  : PredictorBase(conf, jobs) {
}

// -----------------------------------------------------------------------------

void SineOfLatitude::compute(const ioda::ObsSpace & odb,
                             const GeoVaLs &,
                             const ObsDiagnostics &,
                             Eigen::MatrixXd & out) const {
  const std::size_t njobs = jobs_.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.rows() == njobs && out.cols() == nlocs);

  // retrieve the sensor view angle
  std::vector<float> cenlat(nlocs, 0.0);
  odb.get_db("MetaData", "latitude", cenlat);

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    out.col(jl).setConstant(sin(cenlat[jl]*Constants::deg2rad));
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
