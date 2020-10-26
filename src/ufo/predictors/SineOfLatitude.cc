/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "ufo/predictors/SineOfLatitude.h"

#include "ioda/ObsSpace.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<SineOfLatitude>
       makerFuncSineOfLatitude_("sine_of_latitude");

// -----------------------------------------------------------------------------

SineOfLatitude::SineOfLatitude(const eckit::Configuration & conf,
                               const std::vector<int> & jobs,
                               const std::string & sensor,
                               const eckit::mpi::Comm & comm)
  : PredictorBase(conf, jobs, sensor, comm) {
}

// -----------------------------------------------------------------------------

void SineOfLatitude::compute(const ioda::ObsSpace & odb,
                             const GeoVaLs &,
                             const ObsDiagnostics &,
                             ioda::ObsVector & out) const {
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  // retrieve the sensor view angle
  std::vector<float> cenlat(nlocs, 0.0);
  odb.get_db("MetaData", "latitude", cenlat);

  const std::size_t njobs = jobs_.size();
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    for (std::size_t jb = 0; jb < njobs; ++jb) {
      out[jl*njobs+jb] = sin(cenlat[jl] * Constants::deg2rad);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
