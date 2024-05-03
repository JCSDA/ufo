/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "tools/new_obsop/example/ObsExampleTLAD.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/example/ObsExampleTLAD.interface.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsExampleTLAD> makerExampleTL_("Example");
// -----------------------------------------------------------------------------

ObsExampleTLAD::ObsExampleTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_example_tlad_setup_f90(keyOper_, parameters.toConfiguration(), odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsExampleTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsExampleTLAD::~ObsExampleTLAD() {
  ufo_example_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExampleTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & ydiags,
                                   const QCFlags_t & qc_flags) {
  ufo_example_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace(), ydiags.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                   const QCFlags_t & qc_flags) const {
  ufo_example_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                   const QCFlags_t & qc_flags) const {
  ufo_example_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::print(std::ostream & os) const {
  os << "ObsExampleTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
