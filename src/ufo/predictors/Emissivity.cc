/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/predictors/Emissivity.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

#include "oops/base/ObsVariables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

static PredictorMaker<Emissivity> makerFuncEmissivity_("emissivityJacobian");

// -----------------------------------------------------------------------------

Emissivity::Emissivity(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars) {
  // required variables
  geovars_ += oops::Variables({"water_area_fraction"});
  if (vars.size() > 0) {
    hdiags_ += oops::ObsVariables({"brightness_temperature_jacobian_surface_emissivity"},
                               vars.channels());
  } else {
    oops::Log::error() << "Channels size is ZERO !" << std::endl;
    ABORT("Channels size is ZERO !");
  }
}

// -----------------------------------------------------------------------------

void Emissivity::compute(const ioda::ObsSpace & odb,
                         const GeoVaLs & geovals,
                         const ObsDiagnostics & ydiags,
                         const ObsBias &,
                         ioda::ObsVector & out) const {
  const std::size_t nvars = out.nvars();
  const std::size_t nlocs = out.nlocs();

  std::vector<float> pred(nlocs, 0.0);
  std::vector<float> h2o_frac(nlocs, 0.0);
  geovals.get(h2o_frac, oops::Variable{"water_area_fraction"});
  std::string hdiags;
  out.zero();
  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    hdiags = "brightness_temperature_jacobian_surface_emissivity_" +
             std::to_string(vars_.channels()[jvar]);
    ydiags.get(pred, hdiags);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      if (h2o_frac[jloc] < 0.99 && std::fabs(pred[jloc]) > 0.001) {
        out[jloc*nvars+jvar] = pred[jloc];
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
