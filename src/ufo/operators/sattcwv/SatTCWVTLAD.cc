/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/sattcwv/SatTCWVTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<SatTCWVTLAD> makerSatTCWVTL_("SatTCWV");
// -----------------------------------------------------------------------------

SatTCWVTLAD::SatTCWVTLAD(const ioda::ObsSpace & odb,
                         const Parameters_ & params)
  : LinearObsOperatorBase(odb), varin_(), k_matrix(), traj_init(false)
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "surface_pressure"};
  varin_.reset(new oops::Variables(vv));

  oops::Log::trace() << "SatTCWVTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

SatTCWVTLAD::~SatTCWVTLAD() {
  traj_init = false;
  oops::Log::trace() << "SatTCWVTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                const QCFlags_t & qc_flags) {
  // Get number of obs locations & number of model pressure levels
  nprofiles = geovals.nlocs();
  nlevels   = geovals.nlevs(oops::Variable{"air_pressure_levels"});  // number of full (rho) levels

  // (Re)initialise the Jacobian matrix
  if (!k_matrix.empty()) k_matrix.clear();
  k_matrix.resize(nlevels-1, std::vector<double>(nprofiles));

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

  // Check model fields are top-down
  if (plev.front() > plev.back()) {
    // Bottom-up: fail as JEDI should be top-down
    throw eckit::BadValue("model fields must be ordered from the top down", Here());
  }

  // Calculate partial derivatives d(TCWV)/d(q) for each profile element

  // Loop over profiles

  for (size_t prof = 0; prof < nprofiles; ++prof) {
    k_matrix[nlevels - 2][prof] = (ps[prof] - plev[nlevels - 2][prof]) /
                                    Constants::grav;

    // Loop over the rest of the model layers
    for (size_t lev = 0; lev < nlevels - 2; ++lev) {
      k_matrix[lev][prof] = (plev[lev+1][prof] - plev[lev][prof]) /
                            Constants::grav;
    }
  }
  traj_init = true;
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & hofx, const QCFlags_t & qc_flags) const {
  // Ensure trajectory has already been calculated
  ASSERT(traj_init);

  // Check dimensions against trajectory
  ASSERT(geovals.nlocs() == nprofiles);
  ASSERT(geovals.nlevs(oops::Variable{"air_pressure_levels"}) == nlevels);

  // Check hofx size and initialise hofx to zero
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();

  // Get 3-D specific humidity increments on theta levels (kg/kg), one level at a time
  std::vector<std::vector<double>> q_d(nlevels-1, std::vector<double>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels-1; ++lev) {
    geovals.getAtLevel(q_d[lev], oops::Variable{"specific_humidity"}, lev);
  }

  // Loop through the obs, calculating the increment to the observation hofx
  for (size_t prof = 0; prof < nprofiles; ++prof) {
    for (size_t lev = 0; lev < nlevels-1; ++lev) {
      hofx[prof] += k_matrix[lev][prof] * q_d[lev][prof];
    }
  }
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & hofx, const QCFlags_t & qc_flags) const {
  // Ensure trajectory has already been calculated
  ASSERT(traj_init);

  // Check dimensions against trajectory
  ASSERT(geovals.nlocs() == nprofiles);
  ASSERT(geovals.nlevs(oops::Variable{"air_pressure_levels"}) == nlevels);

  // Check hofx size
  ASSERT(geovals.nlocs() == hofx.nlocs());

  // Get 3-D specific humidity increments on theta levels (kg/kg), one level at a time
  std::vector<std::vector<double>> q_d(nlevels-1, std::vector<double>(nprofiles));
  for (std::size_t lev = 0; lev < nlevels-1; ++lev) {
    geovals.getAtLevel(q_d[lev], oops::Variable{"specific_humidity"}, lev);
  }

  // Get the missing value indicator
  const double missing = util::missingValue<double>();

  // Loop through the obs, adding the increment to the model state
  for (size_t prof = 0; prof < nprofiles; ++prof) {
    if (hofx[prof] != missing) {
      for (size_t lev = 0; lev < nlevels-1; ++lev) {
        q_d[lev][prof] += k_matrix[lev][prof] * hofx[prof];
      }
    }
  }

  // Store the updated model state increments
  for (std::size_t lev = 0; lev < nlevels-1; ++lev) {
    geovals.putAtLevel(q_d[lev], oops::Variable{"specific_humidity"}, lev);
  }
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::print(std::ostream & os) const {
  os << "SatTCWVTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
