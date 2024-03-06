/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/vertinterp/ObsVertInterp.h"

#include <ostream>
#include <vector>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/operators/vertinterp/ObsVertInterp.interface.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsVertInterp> makerVertInterp_("VertInterp");
// -----------------------------------------------------------------------------

ObsVertInterp::ObsVertInterp(const ioda::ObsSpace & odb,
                                   const Parameters_ & params)
  : ObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    keyOperVertInterp_(0), odb_(odb), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices);

  varin_ += nameMap_.convertName(operatorVars_);

  ufo_vertinterp_setup_f90(keyOperVertInterp_, params.toConfiguration(),
                              operatorVars_, operatorVarIndices.data(),
                              operatorVarIndices.size(), varin_);

  oops::Log::trace() << "ObsVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsVertInterp::~ObsVertInterp() {
  ufo_vertinterp_delete_f90(keyOperVertInterp_);
  oops::Log::trace() << "ObsVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                   ObsDiagnostics & odiags,
                                   const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsVertInterp::simulateObs entered" << std::endl;

  ufo_vertinterp_simobs_f90(keyOperVertInterp_, gom.toFortran(), odb_,
                               ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsVertInterp::simulateObs exit" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVertInterp::print(std::ostream & os) const {
  os << "ObsVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
