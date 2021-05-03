/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/Ref/ObsGnssroRefTLAD.h"

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
static LinearObsOperatorMaker<ObsGnssroRefTLAD> makerGnssroRefTL_("GnssroRef");
// -----------------------------------------------------------------------------

ObsGnssroRefTLAD::ObsGnssroRefTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : LinearObsOperatorBase(odb), keyOperGnssroRef_(0), varin_()
{
  const eckit::LocalConfiguration obsOptions(config, "obs options");

  ufo_gnssro_ref_tlad_setup_f90(keyOperGnssroRef_, obsOptions);
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroRefTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroRefTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefTLAD::~ObsGnssroRefTLAD() {
  ufo_gnssro_ref_tlad_delete_f90(keyOperGnssroRef_);
  oops::Log::trace() << "ObsGnssroRefTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                     ObsDiagnostics &) {
  ufo_gnssro_ref_tlad_settraj_f90(keyOperGnssroRef_, geovals.toFortran(), obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_gnssro_ref_simobs_tl_f90(keyOperGnssroRef_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_gnssro_ref_simobs_ad_f90(keyOperGnssroRef_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefTLAD::print(std::ostream & os) const {
  os << "ObsGnssroRefTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
