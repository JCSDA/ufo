/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radarradialvelocity/ObsRadarRadialVelocityTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadarRadialVelocityTLAD>
                makerRadarRadialVelocityTL_("RadarRadialVelocity");
// -----------------------------------------------------------------------------

ObsRadarRadialVelocityTLAD::ObsRadarRadialVelocityTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperRadarRadialVelocity_(0), varin_()
{
  ufo_radarradialvelocity_tlad_setup_f90(keyOperRadarRadialVelocity_,
                                         params.toConfiguration(), odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsRadarRadialVelocityTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarRadialVelocityTLAD::~ObsRadarRadialVelocityTLAD() {
  ufo_radarradialvelocity_tlad_delete_f90(keyOperRadarRadialVelocity_);
  oops::Log::trace() << "ObsRadarRadialVelocityTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocityTLAD::setTrajectory(const GeoVaLs & geovals,
                                               ObsDiagnostics &,
                                               const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsRadarRadialVelocityTLAD::setTrajectory entering" << std::endl;

  ufo_radarradialvelocity_tlad_settraj_f90(keyOperRadarRadialVelocity_, geovals.toFortran(),
                                           obsspace());

  oops::Log::trace() << "ObsRadarRadialVelocityTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocityTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                               const QCFlags_t & qc_flags) const {
  ufo_radarradialvelocity_simobs_tl_f90(keyOperRadarRadialVelocity_, geovals.toFortran(),
                                        obsspace(), ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsRadarRadialVelocityTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocityTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                               const QCFlags_t & qc_flags) const {
  ufo_radarradialvelocity_simobs_ad_f90(keyOperRadarRadialVelocity_, geovals.toFortran(),
                                        obsspace(), ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsRadarRadialVelocityTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocityTLAD::print(std::ostream & os) const {
  os << "ObsRadarRadialVelocityTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
