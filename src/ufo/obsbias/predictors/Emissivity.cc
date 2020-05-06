/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/obsbias/predictors/Emissivity.h"

#include "ioda/ObsSpace.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<Emissivity> makerFuncEmissivity_("emissivity");

// -----------------------------------------------------------------------------

Emissivity::Emissivity(const eckit::Configuration & conf)
  : PredictorBase(conf) {
  // required variables
  this->updateGeovars({"water_area_fraction"});
  this->updateHdiagnostics({"brightness_temperature_jacobian_surface_emissivity_CH"});
}

// -----------------------------------------------------------------------------

void Emissivity::compute(const ioda::ObsSpace & odb,
                         const GeoVaLs & geovals,
                         const ObsDiagnostics & ydiags,
                         const std::vector<int> & jobs,
                         Eigen::MatrixXd & out) const {
  const std::size_t njobs = jobs.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.rows() == njobs && out.cols() == nlocs);

  std::vector <float> pred(nlocs, 0.0);
  std::vector<float> h2o_frac(nlocs, 0.0);
  geovals.get(h2o_frac, "water_area_fraction");
  std::string hdiags;
  for (std::size_t jb = 0; jb < njobs; ++jb) {
    hdiags = "brightness_temperature_jacobian_surface_emissivity_" + std::to_string(jobs[jb]);
    ydiags.get(pred, hdiags);
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      if (h2o_frac[jl] < 0.99 && std::fabs(pred[jl]) > 0.001) {
        out(jb, jl) = pred[jl];
      } else {
        out(jb, jl) = 0.0;
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
