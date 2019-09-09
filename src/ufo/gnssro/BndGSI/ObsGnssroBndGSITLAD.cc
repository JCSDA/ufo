/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/BndGSI/ObsGnssroBndGSITLAD.h"

#include <ostream>
#include <string>
#include <vector>


#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroBndGSITLAD> makerGnssroBndGSITL_("GnssroBndGSI");
// -----------------------------------------------------------------------------

ObsGnssroBndGSITLAD::ObsGnssroBndGSITLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOperGnssroBndGSI_(0), varin_(), odb_(odb)
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
  } else {
    vv.push_back("air_pressure_levels");
  }

  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_bndgsi_tlad_setup_f90(keyOperGnssroBndGSI_, &configc);

  oops::Log::info() << "ObsGnssroBndGSITLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBndGSITLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndGSITLAD::~ObsGnssroBndGSITLAD() {
  ufo_gnssro_bndgsi_tlad_delete_f90(keyOperGnssroBndGSI_);
  oops::Log::trace() << "ObsGnssroBndGSITLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSITLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_gnssro_bndgsi_tlad_settraj_f90(keyOperGnssroBndGSI_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSITLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_gnssro_bndgsi_simobs_tl_f90(keyOperGnssroBndGSI_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSITLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_gnssro_bndgsi_simobs_ad_f90(keyOperGnssroBndGSI_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndGSITLAD::print(std::ostream & os) const {
  os << "ObsGnssroBndGSITLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
