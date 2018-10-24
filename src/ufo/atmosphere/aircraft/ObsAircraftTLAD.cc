/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/aircraft/ObsAircraftTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAircraftTLAD> makerAircraftTL_("Aircraft");
// -----------------------------------------------------------------------------

ObsAircraftTLAD::ObsAircraftTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAircraft_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aircraft_tlad_setup_f90(keyOperAircraft_, &configc);
  const std::vector<std::string> vv{"virtual_temperature"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAircraftTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAircraftTLAD::~ObsAircraftTLAD() {
  ufo_aircraft_tlad_delete_f90(keyOperAircraft_);
  oops::Log::trace() << "ObsAircraftTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAircraftTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_aircraft_tlad_settraj_f90(keyOperAircraft_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsAircraftTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_aircraft_simobs_tl_f90(keyOperAircraft_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAircraftTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_aircraft_simobs_ad_f90(keyOperAircraft_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAircraftTLAD::print(std::ostream & os) const {
  os << "ObsAircraftTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
