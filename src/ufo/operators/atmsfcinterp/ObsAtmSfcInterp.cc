/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/atmsfcinterp/ObsAtmSfcInterp.h"

#include <ostream>
#include <vector>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmSfcInterp> makerGSISfcModel_("GSISfcModel");
// -----------------------------------------------------------------------------

ObsAtmSfcInterp::ObsAtmSfcInterp(const ioda::ObsSpace & odb, const Parameters_ & parameters)
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
    keyOperAtmSfcInterp_(0), odb_(odb), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(parameters.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices);

  varin_ += nameMap_.convertName(operatorVars_);

  ufo_atmsfcinterp_setup_f90(keyOperAtmSfcInterp_, parameters.toConfiguration(),
                             operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(),
                             varin_);

  oops::Log::trace() << "ObsAtmSfcInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmSfcInterp::~ObsAtmSfcInterp() {
  ufo_atmsfcinterp_delete_f90(keyOperAtmSfcInterp_);
  oops::Log::trace() << "ObsAtmSfcInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                  ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_atmsfcinterp_simobs_f90(keyOperAtmSfcInterp_, gom.toFortran(), odb_,
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAtmSfcInterp: observation operator executed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::print(std::ostream & os) const {
  os << "ObsAtmSfcInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
