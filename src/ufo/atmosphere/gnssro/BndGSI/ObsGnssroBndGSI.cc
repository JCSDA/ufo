/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/gnssro/BndGSI/ObsGnssroBndGSI.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndGSI> makerGnssroBndGSI_("GnssroBndGSI");
// -----------------------------------------------------------------------------

ObsGnssroBndGSI::ObsGnssroBndGSI(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperGnssroBndGSI_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vv{"temperature", "specific_humidity", "air_pressure",
                                    "air_pressure_levels", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"bending_angle"};
  varout_.reset(new oops::Variables(vout));

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;

  ufo_gnssro_bndgsi_setup_f90(keyOperGnssroBndGSI_, &configc);

  oops::Log::trace() << "ObsGnssroBndGSI created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndGSI::~ObsGnssroBndGSI() {
  ufo_gnssro_bndgsi_delete_f90(keyOperGnssroBndGSI_);
  oops::Log::trace() << "ObsGnssroBndGSI destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSI::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_gnssro_bndgsi_simobs_f90(keyOperGnssroBndGSI_, gom.toFortran(), odb_,
                               ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

Locations * ObsGnssroBndGSI::locateObs(const util::DateTime & t1,
                                       const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_gnssro_bndgsi_locateobs_f90(keyOperGnssroBndGSI_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSI::print(std::ostream & os) const {
  os << "ObsGnssroBndGSI::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
