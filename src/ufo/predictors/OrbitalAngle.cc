/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/predictors/OrbitalAngle.h"
#include "ufo/utils/Constants.h"

namespace ufo {

constexpr char FourierTermTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<FourierTermType>
  FourierTermTypeParameterTraitsHelper::namedValues[];

static PredictorMaker<OrbitalAngle> makerFuncOrbitalAngle_("satelliteOrbitalAngle");

// -----------------------------------------------------------------------------

OrbitalAngle::OrbitalAngle(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order),
    component_(parameters.component) {
  // override the predictor name to distinguish between Orbital angle predictors of
  // different orders as well as the two components, sine and cosine.
  name() = name() + "_order_" + std::to_string(order_)+ "_" +
      (component_ == FourierTermType::SIN ? "sin" : "cos");
}

// -----------------------------------------------------------------------------

void OrbitalAngle::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        const ObsBias &,
                        ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  // retrieve the sensor orbital angle
  std::vector<double> orbital_angle(nlocs, 0.0);
  odb.get_db("MetaData", "satelliteOrbitalAngle", orbital_angle);

  switch (component_)
  {
    case FourierTermType::COS:
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        double cos_oa{ std::cos(orbital_angle[jl]*order_*Constants::deg2rad)};
        for (std::size_t jb = 0; jb < nvars; ++jb) {
          out[jl*nvars+jb] = cos_oa;
        }
      }
      break;
    case FourierTermType::SIN:
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
