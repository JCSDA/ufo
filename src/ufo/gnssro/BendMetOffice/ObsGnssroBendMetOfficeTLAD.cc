/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/BendMetOffice/ObsGnssroBendMetOfficeTLAD.h"

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
static LinearObsOperatorMaker<ObsGnssroBendMetOfficeTLAD>
    makerGnssroBendMetOfficeTL_("GnssroBendMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroBendMetOfficeTLAD::ObsGnssroBendMetOfficeTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOperGnssroBendMetOffice_(0), odb_(odb), varin_()
{
  const eckit::LocalConfiguration obsOptions(config, "obs options");
  const eckit::Configuration * configc = &obsOptions;

  ufo_gnssro_bendmetoffice_tlad_setup_f90(keyOperGnssroBendMetOffice_, &configc);
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroBendMetOfficeTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBendMetOfficeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBendMetOfficeTLAD::~ObsGnssroBendMetOfficeTLAD() {
  ufo_gnssro_bendmetoffice_tlad_delete_f90(keyOperGnssroBendMetOffice_);
  oops::Log::trace() << "ObsGnssroBendMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                           ObsDiagnostics &) {
  ufo_gnssro_bendmetoffice_tlad_settraj_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_gnssro_bendmetoffice_simobs_tl_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_gnssro_bendmetoffice_simobs_ad_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBendMetOfficeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
