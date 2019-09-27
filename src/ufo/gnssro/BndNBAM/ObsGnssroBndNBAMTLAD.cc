/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/BndNBAM/ObsGnssroBndNBAMTLAD.h"

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
static LinearObsOperatorMaker<ObsGnssroBndNBAMTLAD> makerGnssroBndNBAMTL_("GnssroBndNBAM");
// -----------------------------------------------------------------------------

ObsGnssroBndNBAMTLAD::ObsGnssroBndNBAMTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOperGnssroBndNBAM_(0), odb_(odb), varin_()
{
  std::vector<std::string> vv{"air_temperature", "specific_humidity"};

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;

  std::string vertlayer;

//---- get vertical coordinate from config ------------------------
  vertlayer = obsOptions.getString("vertlayer", "full");

  if ( vertlayer == "mass" ) {
    vv.push_back("air_pressure");
  } else {
    vv.push_back("air_pressure_levels");
  }

  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_bndnbam_tlad_setup_f90(keyOperGnssroBndNBAM_, &configc);

  oops::Log::info() << "ObsGnssroBndNBAMTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBndNBAMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndNBAMTLAD::~ObsGnssroBndNBAMTLAD() {
  ufo_gnssro_bndnbam_tlad_delete_f90(keyOperGnssroBndNBAM_);
  oops::Log::trace() << "ObsGnssroBndNBAMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAMTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_gnssro_bndnbam_tlad_settraj_f90(keyOperGnssroBndNBAM_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_gnssro_bndnbam_simobs_tl_f90(keyOperGnssroBndNBAM_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_gnssro_bndnbam_simobs_ad_f90(keyOperGnssroBndNBAM_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAMTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBndNBAMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
