/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/insitutemperature/ObsInsituTemperature.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsInsituTemperature> makerInsituTemperature_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperature::ObsInsituTemperature(const ioda::ObsSpace & odb,
                                           const ObsInsituTemperatureParameters & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  // get optional list of variables to operate on (to be consistent with other
  // operators that work under "composite") eventhough we assume on the Fortran
  // side that this operator will operate on just a single variable (sea_water_temperature)
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);

  // sanity check to make sure waterTemperature is the ONLY variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "waterTemperature",
    "InsituTemperature can only work on variable \"waterTemperature\"");

  ufo_insitutemperature_setup_f90(keyOper_, params.toConfiguration(),
    operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(), varin_);

  oops::Log::trace() << "ObsInsituTemperature created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperature::~ObsInsituTemperature() {
  ufo_insitutemperature_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituTemperature destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                       ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_insitutemperature_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(),
                                   ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperature: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::print(std::ostream & os) const {
  os << "ObsInsituTemperature::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
