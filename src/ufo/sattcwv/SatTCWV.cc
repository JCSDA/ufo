/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/sattcwv/SatTCWV.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// ----------------------------------------------------------------------------
static ObsOperatorMaker<SatTCWV> makerSatTCWV_("SatTCWV");
// -----------------------------------------------------------------------------

SatTCWV::SatTCWV(const ioda::ObsSpace & odb,
                 const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "surface_pressure"};
  varin_.reset(new oops::Variables(vv));

  oops::Log::trace() << "SatTCWV created." << std::endl;
}

// -----------------------------------------------------------------------------

SatTCWV::~SatTCWV() {
  oops::Log::trace() << "SatTCWV destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SatTCWV::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                          ObsDiagnostics &) const {
  // Check hofx size and initialise hofx to zero
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();

  // Get number of obs locations & number of model pressure levels
  std::size_t nprofiles = geovals.nlocs();
  std::size_t nlevels   = geovals.nlevs("air_pressure_levels");  // number of full (rho) levels

  // Get 2-D surface pressure
  std::vector<float> ps(nprofiles);  // surface pressure (Pa)
  geovals.get(ps, "surface_pressure");

  // Get 3-D air pressure on rho levels (Pa), one level at a time
  std::vector<std::vector<float>> plev(nlevels, std::vector<float>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels; ++lev) {
    geovals.getAtLevel(plev[lev], "air_pressure_levels", lev);
  }

  // Get 3-D specific humidity on theta levels (kg/kg), one level at a time
  std::vector<std::vector<float>> q(nlevels-1, std::vector<float>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels-1; ++lev) {
    geovals.getAtLevel(q[lev], "specific_humidity", lev);
  }

  // Check whether model fields are top-down or bottom-up
  int lev_near_surf, lay_near_surf, lev1, lev2;
  if (plev[2][0] - plev[1][0] > 0) {
    // Top-down
    lev_near_surf = nlevels - 2;
    lay_near_surf = nlevels - 2;
    lev1 = 0;
    lev2 = nlevels - 2;
  } else {
    // Bottom-up
    lev_near_surf = 1;
    lay_near_surf = 0;
    lev1 = 1;
    lev2 = nlevels - 1;
  }

  // Calculate TCWV for each profile, integrating over each layer

  // Loop over profiles
  for (size_t prof = 0; prof < nprofiles; ++prof) {
    // Start with lowest model layer, using surface pressure
    // NB this assumes surface q is same as q 10m but could use q2m in future
    hofx[prof] = abs(ps[prof] - plev[lev_near_surf][prof]) *
                 q[lay_near_surf][prof] / Constants::grav;

    // Loop over the rest of the model layers
    for (size_t lev = lev1; lev < lev2; ++lev) {
      hofx[prof] += abs(plev[lev][prof] - plev[lev+1][prof]) *
                    q[lev][prof] / Constants::grav;
    }
  }
}

// -----------------------------------------------------------------------------

void SatTCWV::print(std::ostream & os) const {
  os << "SatTCWV::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
