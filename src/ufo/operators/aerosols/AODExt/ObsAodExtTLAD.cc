/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/AODExt/ObsAodExtTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodExtTLAD> makerAodExtTL_("AodExt");
// -----------------------------------------------------------------------------

ObsAodExtTLAD::ObsAodExtTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_aodext_tlad_setup_f90(keyOper_, params.toConfiguration(),
                            odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsAodExtTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodExtTLAD::~ObsAodExtTLAD() {
  ufo_aodext_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodExtTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExtTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                  const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsAodExtTLAD:trajectory entering" <<std::endl;
  ufo_aodext_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsAodExtTLAD: set trajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExtTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                  const QCFlags_t & qc_flags) const {
  ufo_aodext_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodExtTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExtTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                  const QCFlags_t & qc_flags) const {
  ufo_aodext_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodExtTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExtTLAD::print(std::ostream & os) const {
  os << "ObsAodExtTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
