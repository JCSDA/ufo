/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/BndROPP2D/ObsGnssroBndROPP2DTLAD.h"

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
static LinearObsOperatorMaker<ObsGnssroBndROPP2DTLAD> makerGnssroBndROPP2DTL_("GnssroBndROPP2D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP2DTLAD::ObsGnssroBndROPP2DTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOperGnssroBndROPP2D_(0), varin_(), odb_(odb)
{
  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;

  ufo_gnssro_bndropp2d_tlad_setup_f90(keyOperGnssroBndROPP2D_, &configc);
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroBndROPP2DTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBndROPP2DTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP2DTLAD::~ObsGnssroBndROPP2DTLAD() {
  ufo_gnssro_bndropp2d_tlad_delete_f90(keyOperGnssroBndROPP2D_);
  oops::Log::trace() << "ObsGnssroBndROPP2DTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_gnssro_bndropp2d_tlad_settraj_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_gnssro_bndropp2d_simobs_tl_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_gnssro_bndropp2d_simobs_ad_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP2DTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
