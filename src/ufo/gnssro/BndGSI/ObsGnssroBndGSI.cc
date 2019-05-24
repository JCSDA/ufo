/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/gnssro/BndGSI/ObsGnssroBndGSI.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndGSI> makerGnssroBndGSI_("GnssroBndGSI");
// -----------------------------------------------------------------------------

ObsGnssroBndGSI::ObsGnssroBndGSI(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperGnssroBndGSI_(0), odb_(odb), varin_(), varout_()
{
  std::vector<std::string> vv{"air_temperature", "specific_humidity"};

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;

  std::string vertlayer;

//---- get vertical coordinate from config ------------------------
  if ( obsOptions.has("vertlayer") ) {
     vertlayer = obsOptions.getString("vertlayer");
  } else {
     vertlayer = "full";
  }

  if ( vertlayer == "mass" ) {
    vv.push_back("air_pressure");
    vv.push_back("geopotential_height");
  } else {
    vv.push_back("air_pressure_levels");
    vv.push_back("geopotential_height_levels");
  }

  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"bending_angle"};
  varout_.reset(new oops::Variables(vout));

  ufo_gnssro_bndgsi_setup_f90(keyOperGnssroBndGSI_, &configc);
  oops::Log::trace() << "ObsGnssroBndGSI created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndGSI::~ObsGnssroBndGSI() {
  ufo_gnssro_bndgsi_delete_f90(keyOperGnssroBndGSI_);
  oops::Log::trace() << "ObsGnssroBndGSI destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSI::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec) const {
  ufo_gnssro_bndgsi_simobs_f90(keyOperGnssroBndGSI_, gom.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSI::print(std::ostream & os) const {
  os << "ObsGnssroBndGSI::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
