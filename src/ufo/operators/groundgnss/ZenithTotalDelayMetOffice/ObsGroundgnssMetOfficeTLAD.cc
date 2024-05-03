/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/groundgnss/ZenithTotalDelayMetOffice/ObsGroundgnssMetOfficeTLAD.h"

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
static LinearObsOperatorMaker<ObsGroundgnssMetOfficeTLAD>
    makerGroundgnssMetOfficeTL_("GroundgnssMetOffice");
// -----------------------------------------------------------------------------

ObsGroundgnssMetOfficeTLAD::ObsGroundgnssMetOfficeTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOperGroundgnssMetOffice_(0), varin_()
{
  ufo_groundgnss_metoffice_tlad_setup_f90(keyOperGroundgnssMetOffice_,
                                          parameters.toConfiguration());
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGroundgnssMetOfficeTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssMetOfficeTLAD::~ObsGroundgnssMetOfficeTLAD() {
  ufo_groundgnss_metoffice_tlad_delete_f90(keyOperGroundgnssMetOffice_);
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                               const QCFlags_t & qc_flags) {
  ufo_groundgnss_metoffice_tlad_settraj_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                            obsspace());
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_groundgnss_metoffice_simobs_tl_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_groundgnss_metoffice_simobs_ad_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsGroundgnssMetOfficeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
