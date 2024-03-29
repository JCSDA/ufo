/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/marinevertinterp/ObsMarineVertInterpTLAD.h"

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsMarineVertInterpTLAD> makerMarinevertinterpTL_("MarineVertInterp");
// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::ObsMarineVertInterpTLAD(const ioda::ObsSpace & odb,
                                                 const ObsMarineVertInterpParameters & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);

  ufo_marinevertinterp_tlad_setup_f90(keyOper_, params.toConfiguration(),
    operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(), varin_);

  oops::Log::trace() << "ObsMarineVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::~ObsMarineVertInterpTLAD() {
  ufo_marinevertinterp_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &) {
  ufo_marinevertinterp_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_marinevertinterp_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                                     ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_marinevertinterp_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                                     ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsMarinevertinterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
