/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/predictors/OrbitalAngle.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<OrbitalAngle> makerFuncOrbitalAngle_("orbital_angle");

// -----------------------------------------------------------------------------

OrbitalAngle::OrbitalAngle(const eckit::Configuration & conf, const oops::Variables & vars)
  : PredictorBase(conf, vars), order_(1) {
    // get the order if it is provided in options
    conf.get("options.order", order_);
    conf.get("options.component", component_);

    // override the predictor name to distinguish between Orbital angle predictors of
    // different orders as well as the two components, sine and cosine.
    name() = name() + "_" + std::to_string(order_)+ "_" + component_;
}

// -----------------------------------------------------------------------------

void OrbitalAngle::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  // retrieve the sensor orbital angle
  std::vector<double> orbital_angle(nlocs, 0.0);
  odb.get_db("MetaData", "satellite_orbital_angle", orbital_angle);

  ASSERT(component_ == "cos" || component_ == "sin");
  switch (component_ == "cos")
  {
    case true:
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        double cos_oa{ std::cos(orbital_angle[jl]*order_*Constants::deg2rad)};
        for (std::size_t jb = 0; jb < nvars; ++jb) {
          out[jl*nvars+jb] = cos_oa;
        }
      }
      break;
    case false:
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        double sin_oa{ std::sin(orbital_angle[jl]*order_*Constants::deg2rad)};
        for (std::size_t jb = 0; jb < nvars; ++jb) {
          out[jl*nvars+jb] = sin_oa;
        }
      }
      break;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
