/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOfficeParameters.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOfficeTLAD.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroRefMetOfficeTLAD>
    makerGnssroRefMetOfficeTL_("GnssroRefMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroRefMetOfficeTLAD::ObsGnssroRefMetOfficeTLAD(const ioda::ObsSpace & odb,
                                                     const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOperGnssroRefMetOffice_(0), varin_(),
    parameters_(parameters)
{
  ObsGnssroRefMetOfficeOptions obsOptions = parameters_.obsOptions.value();

  ufo_gnssro_refmetoffice_tlad_setup_f90(keyOperGnssroRefMetOffice_,
                                         obsOptions.vertInterpOPS,
                                         obsOptions.pseudoLevels,
                                         obsOptions.minTempGrad);
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroRefMetOfficeTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroRefMetOfficeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefMetOfficeTLAD::~ObsGnssroRefMetOfficeTLAD() {
  ufo_gnssro_refmetoffice_tlad_delete_f90(keyOperGnssroRefMetOffice_);
  oops::Log::trace() << "ObsGnssroRefMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                              const QCFlags_t & qc_flags) {
  ufo_gnssro_refmetoffice_tlad_settraj_f90(keyOperGnssroRefMetOffice_, geovals.toFortran(),
                                           obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOfficeTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_gnssro_refmetoffice_simobs_tl_f90(keyOperGnssroRefMetOffice_, geovals.toFortran(),
                                        obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOfficeTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_gnssro_refmetoffice_simobs_ad_f90(keyOperGnssroRefMetOffice_, geovals.toFortran(),
                                        obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsGnssroRefMetOfficeTLAD: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
