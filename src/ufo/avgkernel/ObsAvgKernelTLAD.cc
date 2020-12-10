/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/avgkernel/ObsAvgKernelTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAvgKernelTLAD> makerAvgKernelTL_("AvgKernel");
// -----------------------------------------------------------------------------

ObsAvgKernelTLAD::ObsAvgKernelTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOperAvgKernel_(0), odb_(odb), varin_()
{
  ufo_avgkernel_tlad_setup_f90(keyOperAvgKernel_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsAvgKernelTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAvgKernelTLAD::~ObsAvgKernelTLAD() {
  ufo_avgkernel_tlad_delete_f90(keyOperAvgKernel_);
  oops::Log::trace() << "ObsAvgKernelTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernelTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                         ObsDiagnostics &) {
  oops::Log::trace() << "ObsAvgKernelTLAD::setTrajectory entering" << std::endl;

  ufo_avgkernel_tlad_settraj_f90(keyOperAvgKernel_, geovals.toFortran(), odb_);

  oops::Log::trace() << "ObsAvgKernelTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernelTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_avgkernel_simobs_tl_f90(keyOperAvgKernel_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAvgKernelTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernelTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_avgkernel_simobs_ad_f90(keyOperAvgKernel_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAvgKernelTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernelTLAD::print(std::ostream & os) const {
  os << "ObsAvgKernelTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
