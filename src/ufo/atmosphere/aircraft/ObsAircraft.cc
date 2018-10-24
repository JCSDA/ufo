/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/aircraft/ObsAircraft.h"

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
static ObsOperatorMaker<ObsAircraft> makerAircraft_("Aircraft");
// -----------------------------------------------------------------------------

ObsAircraft::ObsAircraft(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAircraft_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  const eckit::Configuration * configc = &config;
  ufo_aircraft_setup_f90(keyOperAircraft_, &configc);
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAircraft created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAircraft::~ObsAircraft() {
  ufo_aircraft_delete_f90(keyOperAircraft_);
  oops::Log::trace() << "ObsAircraft destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAircraft::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_aircraft_simobs_f90(keyOperAircraft_, gom.toFortran(), odb_,
                          ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAircraft::print(std::ostream & os) const {
  os << "ObsAircraft::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
