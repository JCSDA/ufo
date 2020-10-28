/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "ufo/predictors/Emissivity.h"

#include "ioda/ObsSpace.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

static PredictorMaker<Emissivity> makerFuncEmissivity_("emissivity");

// -----------------------------------------------------------------------------

Emissivity::Emissivity(const eckit::Configuration & conf, const std::vector<int> & jobs)
  : PredictorBase(conf, jobs) {
  // required variables
  geovars_ += oops::Variables({"water_area_fraction"});
  if (jobs.size() > 0) {
    hdiags_ += oops::Variables({"brightness_temperature_jacobian_surface_emissivity"}, jobs);
  } else {
    oops::Log::error() << "Channels size is ZERO !" << std::endl;
    ABORT("Channels size is ZERO !");
  }
}

// -----------------------------------------------------------------------------

void Emissivity::compute(const ioda::ObsSpace & odb,
                         const GeoVaLs & geovals,
                         const ObsDiagnostics & ydiags,
                         ioda::ObsVector & out) const {
  const std::size_t njobs = jobs_.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  std::vector <float> pred(nlocs, 0.0);
  std::vector<float> h2o_frac(nlocs, 0.0);
  geovals.get(h2o_frac, "water_area_fraction");
  std::string hdiags;
  out.zero();
  for (std::size_t jb = 0; jb < njobs; ++jb) {
    hdiags = "brightness_temperature_jacobian_surface_emissivity_" + std::to_string(jobs_[jb]);
    ydiags.get(pred, hdiags);
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      if (h2o_frac[jl] < 0.99 && std::fabs(pred[jl]) > 0.001) {
        out[jl*njobs+jb] = pred[jl];
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
