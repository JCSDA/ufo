/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/identity/ObsIdentityTLAD.h"

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsIdentityTLAD> makerIdentityTL_("Identity");
// -----------------------------------------------------------------------------

ObsIdentityTLAD::ObsIdentityTLAD(const ioda::ObsSpace & odb,
                                 const eckit::Configuration & config)
  : LinearObsOperatorBase(odb), keyOperObsIdentity_(0), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(config, odb.obsvariables(), operatorVars_, operatorVarIndices);

  ufo_identity_tlad_setup_f90(keyOperObsIdentity_, config,
                              operatorVars_,
                              operatorVarIndices.data(), operatorVarIndices.size(),
                              varin_);

  oops::Log::trace() << "ObsIdentityTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsIdentityTLAD::~ObsIdentityTLAD() {
  ufo_identity_tlad_delete_f90(keyOperObsIdentity_);
  oops::Log::trace() << "ObsIdentityTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentityTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                    ObsDiagnostics &) {
  ufo_identity_tlad_settraj_f90(keyOperObsIdentity_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsIdentityTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentityTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_identity_simobs_tl_f90(keyOperObsIdentity_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsIdentityTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentityTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_identity_simobs_ad_f90(keyOperObsIdentity_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsIdentityTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentityTLAD::print(std::ostream & os) const {
  os << "ObsIdentityTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
