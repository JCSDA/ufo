/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/marinevertinterp/ObsMarineVertInterp.h"

#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables
namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsMarineVertInterp> makerMarineVertInterp_("MarineVertInterp");
// -----------------------------------------------------------------------------

ObsMarineVertInterp::ObsMarineVertInterp(const ioda::ObsSpace & odb,
                                         const ObsMarineVertInterpParameters & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);

  ufo_marinevertinterp_setup_f90(keyOper_, params.toConfiguration(),
    operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(),  varin_);

  oops::Log::trace() << "ObsMarineVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterp::~ObsMarineVertInterp() {
  ufo_marinevertinterp_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                      ObsDiagnostics &) const {
  ufo_marinevertinterp_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(),
                                  ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::print(std::ostream & os) const {
  os << "ObsMarineVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
