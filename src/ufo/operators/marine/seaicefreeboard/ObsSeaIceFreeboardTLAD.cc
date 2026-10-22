/*
 * (C) Copyright 2026 NOAA/NWS/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/seaicefreeboard/ObsSeaIceFreeboardTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceFreeboardTLAD> makerSeaIceFreeboardTL_("SeaIceFreeboard");
// -----------------------------------------------------------------------------

ObsSeaIceFreeboardTLAD::ObsSeaIceFreeboardTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), varin_(seaicefreeboard::requiredVars()),
    coefs_(seaicefreeboard::coefficients(params))
{
  seaicefreeboard::checkSimulatedVariable(odb);
  oops::Log::trace() << "ObsSeaIceFreeboardTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFreeboardTLAD::~ObsSeaIceFreeboardTLAD() {
  oops::Log::trace() << "ObsSeaIceFreeboardTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboardTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                           const QCFlags_t &) {
  const double missing = util::missingValue<double>();
  const size_t nlocs = geovals.nlocs();

  seaicefreeboard::checkSingleCategory(geovals);

  // The operator is linear in the thicknesses, so only the missing-value pattern
  // of the trajectory is needed, not its values.
  std::vector<double> iceThick(nlocs);
  std::vector<double> snowThick(nlocs);
  geovals.getAtLevel(iceThick, oops::Variable{"sea_ice_thickness"}, 0);
  geovals.getAtLevel(snowThick, oops::Variable{"sea_ice_snow_thickness"}, 0);

  missingTraj_.assign(nlocs, false);
  for (size_t i = 0; i < nlocs; ++i) {
    missingTraj_[i] = (iceThick[i] == missing || snowThick[i] == missing);
  }

  oops::Log::trace() << "ObsSeaIceFreeboardTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboardTLAD::simulateObsTL(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  const double missing = util::missingValue<double>();
  const size_t nlocs = ovec.size();

  std::vector<double> dIceThick(nlocs);
  std::vector<double> dSnowThick(nlocs);
  gv.getAtLevel(dIceThick, oops::Variable{"sea_ice_thickness"}, 0);
  gv.getAtLevel(dSnowThick, oops::Variable{"sea_ice_snow_thickness"}, 0);

  for (size_t i = 0; i < nlocs; ++i) {
    if (missingTraj_[i]) {
      ovec[i] = missing;
    } else {
      ovec[i] = coefs_.alpha * dIceThick[i] + coefs_.beta * dSnowThick[i];
    }
  }

  oops::Log::trace() << "ObsSeaIceFreeboardTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboardTLAD::simulateObsAD(GeoVaLs & gv, const ioda::ObsVector & ovec) const {
  const double missing = util::missingValue<double>();
  const size_t nlocs = ovec.size();

  // Accumulate into the current adjoint values, which may hold contributions
  // from other operators sharing these geovals.
  std::vector<double> dIceThick(nlocs);
  std::vector<double> dSnowThick(nlocs);
  gv.getAtLevel(dIceThick, oops::Variable{"sea_ice_thickness"}, 0);
  gv.getAtLevel(dSnowThick, oops::Variable{"sea_ice_snow_thickness"}, 0);

  for (size_t i = 0; i < nlocs; ++i) {
    if (ovec[i] != missing && !missingTraj_[i]) {
      dIceThick[i] += coefs_.alpha * ovec[i];
      dSnowThick[i] += coefs_.beta * ovec[i];
    }
  }
  gv.putAtLevel(dIceThick, oops::Variable{"sea_ice_thickness"}, 0);
  gv.putAtLevel(dSnowThick, oops::Variable{"sea_ice_snow_thickness"}, 0);

  oops::Log::trace() << "ObsSeaIceFreeboardTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboardTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFreeboardTLAD operator, alpha=" << coefs_.alpha << ", beta=" << coefs_.beta;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
