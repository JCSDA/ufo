/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/sattcwv/SatTCWV.h"

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
                 const Parameters_ & params)
        : ObsOperatorBase(odb), varin_()
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
                          ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  // Check hofx size and initialise hofx to zero
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();

  // Get number of obs locations & number of model pressure levels
  std::size_t nprofiles = geovals.nlocs();
  // number of full (rho) levels
  std::size_t nlevels   = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  // Get 2-D surface pressure
  std::vector<float> ps(nprofiles);  // surface pressure (Pa)
  geovals.get(ps, oops::Variable{"surface_pressure"});

  // Get 3-D air pressure on rho levels (Pa), one level at a time
  std::vector<std::vector<float>> plev(nlevels, std::vector<float>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels; ++lev) {
    geovals.getAtLevel(plev[lev], oops::Variable{"air_pressure_levels"}, lev);
  }

  // Get 3-D specific humidity on theta levels (kg/kg), one level at a time
  std::vector<std::vector<float>> q(nlevels-1, std::vector<float>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels-1; ++lev) {
    geovals.getAtLevel(q[lev], oops::Variable{"specific_humidity"}, lev);
  }

  // Check model fields are top-down, fail if not
  if (plev.front() > plev.back()) {
    throw eckit::BadValue("model fields must be ordered from the top down", Here());
  }

  // Calculate TCWV for each profile, integrating over each layer

  // Loop over profiles
  for (size_t prof = 0; prof < nprofiles; ++prof) {
    // Start with lowest model layer, using surface pressure
    // NB this assumes surface q is same as q 10m but could use q2m in future
    hofx[prof] = (ps[prof] - plev[nlevels - 2][prof]) *
                 q[nlevels - 2][prof] / Constants::grav;

    // Loop over the rest of the model layers
    for (size_t lev = 0; lev < nlevels - 2; ++lev) {
      hofx[prof] += (plev[lev+1][prof] - plev[lev][prof]) *
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
