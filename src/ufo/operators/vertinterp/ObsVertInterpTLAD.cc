/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/vertinterp/ObsVertInterpTLAD.h"

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/vertinterp/ObsVertInterpTLAD.interface.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsVertInterpTLAD> makerVertInterpTL_("VertInterp");
// -----------------------------------------------------------------------------

ObsVertInterpTLAD::ObsVertInterpTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    keyOperVertInterp_(0), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices);

  varin_ += nameMap_.convertName(operatorVars_);

  ufo_vertinterp_tlad_setup_f90(keyOperVertInterp_, params.toConfiguration(),
                                   operatorVars_,
                                   operatorVarIndices.data(), operatorVarIndices.size(),
                                   varin_);

  oops::Log::trace() << "ObsVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVertInterpTLAD::~ObsVertInterpTLAD() {
  ufo_vertinterp_tlad_delete_f90(keyOperVertInterp_);
  oops::Log::trace() << "ObsVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                      const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsVertInterpTLAD::setTrajectory entering" << std::endl;

  ufo_vertinterp_tlad_settraj_f90(keyOperVertInterp_, geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsVertInterpTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                      const QCFlags_t & qc_flags) const {
  ufo_vertinterp_simobs_tl_f90(keyOperVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsVertInterpTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                      const QCFlags_t & qc_flags) const {
  ufo_vertinterp_simobs_ad_f90(keyOperVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsVertInterpTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsVertInterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
