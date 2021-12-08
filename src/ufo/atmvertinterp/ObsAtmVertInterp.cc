/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmvertinterp/ObsAtmVertInterp.h"

#include <ostream>
#include <vector>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/atmvertinterp/ObsAtmVertInterp.interface.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmVertInterp> makerVertInterp_("VertInterp");
// -----------------------------------------------------------------------------

ObsAtmVertInterp::ObsAtmVertInterp(const ioda::ObsSpace & odb,
                                   const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperAtmVertInterp_(0),
    odb_(odb), varin_()
{
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.obsvariables(),
                       operatorVars_, operatorVarIndices);

  ufo_atmvertinterp_setup_f90(keyOperAtmVertInterp_, params.toConfiguration(),
                              operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(),
                              varin_);

  oops::Log::trace() << "ObsAtmVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterp::~ObsAtmVertInterp() {
  ufo_atmvertinterp_delete_f90(keyOperAtmVertInterp_);
  oops::Log::trace() << "ObsAtmVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                   ObsDiagnostics &) const {
  oops::Log::trace() << "ObsAtmVertInterp::simulateObs entered" << std::endl;

  ufo_atmvertinterp_simobs_f90(keyOperAtmVertInterp_, gom.toFortran(), odb_,
                               ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterp::simulateObs exit" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterp::print(std::ostream & os) const {
  os << "ObsAtmVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
