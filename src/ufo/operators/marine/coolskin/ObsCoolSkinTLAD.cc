/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/coolskin/ObsCoolSkinTLAD.h"

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
static LinearObsOperatorMaker<ObsCoolSkinTLAD> makerCoolSkinTL_("CoolSkin");
// -----------------------------------------------------------------------------

ObsCoolSkinTLAD::ObsCoolSkinTLAD(const ioda::ObsSpace & odb, const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  const std::vector<std::string> vv{"sea_surface_temperature",
                                    "net_downwelling_shortwave_radiation",
                                    "upward_latent_heat_flux_in_air",
                                    "upward_sensible_heat_flux_in_air",
                                    "net_downwelling_longwave_radiation",
                                    "friction_velocity_over_water"};
  varin_.reset(new oops::Variables(vv));
  ufo_CoolSkin_tlad_setup_f90(keyOper_, params.toConfiguration());
  oops::Log::trace() << "ObsCoolSkinTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsCoolSkinTLAD::~ObsCoolSkinTLAD() {
  ufo_CoolSkin_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsCoolSkinTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkinTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                    const QCFlags_t & qc_flags) {
  ufo_CoolSkin_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsCoolSkinTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkinTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
 const QCFlags_t& qc_flags) const {
  ufo_CoolSkin_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                             ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsCoolSkinTLAD: tangent linear observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkinTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    const QCFlags_t & qc_flags) const {
  ufo_CoolSkin_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                             ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsCoolSkinTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkinTLAD::print(std::ostream & os) const {
  os << "ObsCoolSkinTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
