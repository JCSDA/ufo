/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/lightning/ObsLightningTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsLightningTLAD> makerLightningTL_("Lightning");
// -----------------------------------------------------------------------------

ObsLightningTLAD::ObsLightningTLAD(const ioda::ObsSpace & odb, const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  // Set local variables from YAML params.
  nhoriz_ = params.num_gridpoints.value()*params.num_gridpoints.value();
  // l_fed_nonlinear_ = params.use_nonlinear.value();

  // ufo_lightning_tlad_setup_f90(keyOper_, l_fed_nonlinear_, nhoriz_,
  ufo_lightning_tlad_setup_f90(keyOper_, nhoriz_,
                                        odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsLightningTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsLightningTLAD::~ObsLightningTLAD() {
  ufo_lightning_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsLightningTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLightningTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &, const QCFlags_t&) {
  oops::Log::trace() << "ObsLightningTLAD: setTrajectory entering" << std::endl;
  ufo_lightning_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsLightningTLAD: setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLightningTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags) const {
  ufo_lightning_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsLightningTLAD: simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLightningTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags) const {
  ufo_lightning_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                              ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsLightningTLAD: simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLightningTLAD::print(std::ostream & os) const {
  os << "ObsLightningTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
