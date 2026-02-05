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
  const std::vector<std::string> vv{"air_pressure_levels",
      "water_vapor_mixing_ratio_wrt_moist_air"};

  varin_.reset(new oops::Variables(vv));

  ufo_groundgnss_metoffice_tlad_setup_f90(keyOperGroundgnssMetOffice_,
                                          parameters.vertInterpOPS.value(),
                                          parameters.pseudoLevels.value(),
                                          parameters.minTempGrad.value(),
                                          parameters.dryRefractivityConstant.value(),
                                          parameters.wetRefractivityConstant.value());

  oops::Log::info() << "ObsGroundgnssMetOfficeTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssMetOfficeTLAD::~ObsGroundgnssMetOfficeTLAD() {
  ufo_groundgnss_metoffice_tlad_delete_f90(keyOperGroundgnssMetOffice_);
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                               const QCFlags_t & qc_flags) {
  ufo_groundgnss_metoffice_tlad_settraj_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                            obsspace());
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_groundgnss_metoffice_simobs_tl_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD::simulateObsTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_groundgnss_metoffice_simobs_ad_f90(keyOperGroundgnssMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsGroundgnssMetOfficeTLAD::simulateObsAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsGroundgnssMetOfficeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
