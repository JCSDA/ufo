/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/atmvertinterp/ObsAtmVertInterpTLAD.h"

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/atmvertinterp/ObsAtmVertInterpTLAD.interface.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerVertInterpTL_("VertInterp");
// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::ObsAtmVertInterpTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperAtmVertInterp_(0), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices);

  ufo_atmvertinterp_tlad_setup_f90(keyOperAtmVertInterp_, params.toConfiguration(),
                                   operatorVars_,
                                   operatorVarIndices.data(), operatorVarIndices.size(),
                                   varin_);

  oops::Log::trace() << "ObsAtmVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::~ObsAtmVertInterpTLAD() {
  ufo_atmvertinterp_tlad_delete_f90(keyOperAtmVertInterp_);
  oops::Log::trace() << "ObsAtmVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &) {
  oops::Log::trace() << "ObsAtmVertInterpTLAD::setTrajectory entering" << std::endl;

  ufo_atmvertinterp_tlad_settraj_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_atmvertinterp_simobs_tl_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_atmvertinterp_simobs_ad_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsAtmVertInterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
