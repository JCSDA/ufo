/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/LiquidWaterPathGuess.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {
static ObsFunctionMaker<LiquidWaterPathGuess> makerLiquidWaterPathGuess_("LiquidWaterPathGuess");

LiquidWaterPathGuess::LiquidWaterPathGuess(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/mass_content_of_cloud_liquid_water_in_atmosphere_layer");
  invars_ += Variable("GeoVaLs/air_pressure");
}

// -----------------------------------------------------------------------------

LiquidWaterPathGuess::~LiquidWaterPathGuess() {}

// -----------------------------------------------------------------------------

void LiquidWaterPathGuess::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/air_pressure"));

  // define scalar to multiply at end
  const float scalevals = 1.0/(2.0*Constants::grav);


  std::vector<float> lwp(nlocs, 0.0);
  float dp;  // pressure difference across layer (Pa)
  std::vector<float> p_leveli(nlocs);  // pressure at level i  (Pa)
  std::vector<float> p_levelipone(nlocs);  // pressure at level i+1 (Pa)
  std::vector<float> clw_leveli(nlocs);  // clw at level i (kg/kg)
  std::vector<float> clw_levelipone(nlocs);  // clw at level i+1 (kg/kg)

  // initialise  p, clw for top of first layer
  in.get(Variable("GeoVaLs/air_pressure"), 0, p_leveli);
  in.get(Variable("GeoVaLs/mass_content_of_cloud_liquid_water_in_atmosphere_layer"),
    0, clw_leveli);

  // perform LWP (kg/m^2) calculation
  // = sum over layers ( (pi+1-pi) * 0.5* (clwi+clwi+1) / g)
  for (size_t ilev = 0; ilev < nlevs-1; ++ilev) {
    in.get(Variable("GeoVaLs/air_pressure"), ilev+1, p_levelipone);
    in.get(Variable("GeoVaLs/mass_content_of_cloud_liquid_water_in_atmosphere_layer"),
      ilev+1, clw_levelipone);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      dp = std::abs(p_levelipone[iloc] - p_leveli[iloc]);
      lwp[iloc] += (dp * (clw_leveli[iloc]+clw_levelipone[iloc]));
    }
    // copy to next layer in profile. saves extracting from geovals
    p_leveli = p_levelipone;
    clw_leveli = clw_levelipone;
  }
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    lwp[iloc] *= scalevals;
  }
  out[0] = lwp;
}

// -----------------------------------------------------------------------------

const ufo::Variables & LiquidWaterPathGuess::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
