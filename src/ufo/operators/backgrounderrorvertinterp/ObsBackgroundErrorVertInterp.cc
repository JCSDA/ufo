/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/backgrounderrorvertinterp/ObsBackgroundErrorVertInterp.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/backgrounderrorvertinterp/ObsBackgroundErrorVertInterp.interface.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

static ObsOperatorMaker<ObsBackgroundErrorVertInterp> maker("BackgroundErrorVertInterp");

ObsBackgroundErrorVertInterp::ObsBackgroundErrorVertInterp(const ioda::ObsSpace & odb,
                                                           const Parameters_ & parameters)
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
    odb_(odb), parameters_(parameters)
{
  oops::Log::trace() << "ObsBackgroundErrorVertInterp constructor entered" << std::endl;

  requiredVars_.push_back(parameters_.verticalCoordinate);

  /// All simulated variables.
  const oops::ObsVariables & obsVars = odb.assimvariables();

  // If the `variables` option is specified, only the variables in that list will have
  // their background errors computed. Otherwise the background errors for all simulated
  // variables are calculated.
  std::vector<int> operatorVarIndices;
  oops::ObsVariables operatorVars;
  getOperatorVariables(parameters.variables.value(), obsVars,
                       operatorVars, operatorVarIndices);
  for (auto ivar : operatorVarIndices)
    requiredVars_.push_back(nameMap_.convertName(obsVars[ivar]).name() + "_background_error");

  oops::Log::trace() << "ObsBackgroundErrorVertInterp created" << std::endl;
}

ObsBackgroundErrorVertInterp::~ObsBackgroundErrorVertInterp() {
  oops::Log::trace() << "ObsBackgroundErrorVertInterp destructed" << std::endl;
}

void ObsBackgroundErrorVertInterp::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                               ObsDiagnostics & ydiags,
                                               const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsBackgroundErrorVertInterp: simulateObs entered" << std::endl;

  const std::string &obsVerticalCoordinate = parameters_.observationVerticalCoordinate;
  const std::string &obsVerticalGroup = parameters_.observationVerticalGroup;
  const std::string &verticalCoordinate = parameters_.verticalCoordinate;

  oops::Variables variables;
  if (parameters_.variables.value() != boost::none)
    for (const Variable &variable : *parameters_.variables.value())
      variables += nameMap_.convertName(variable.toOopsObsVariables());
  else
    variables = nameMap_.convertName(odb_.assimvariables());

  ufo_backgrounderrorvertinterp_fillobsdiags_f90(obsVerticalCoordinate.size(),
                                                 obsVerticalCoordinate.c_str(),
                                                 obsVerticalGroup.size(),
                                                 obsVerticalGroup.c_str(),
                                                 verticalCoordinate.size(),
                                                 verticalCoordinate.c_str(),
                                                 geovals.toFortran(), odb_, hofx.nlocs(),
                                                 variables,
                                                 ydiags.toFortran());

  oops::Log::trace() << "ObsBackgroundErrorVertInterp: simulateObs exit" <<  std::endl;
}

const oops::Variables & ObsBackgroundErrorVertInterp::requiredVars() const {
  return requiredVars_;
}

oops::ObsVariables ObsBackgroundErrorVertInterp::simulatedVars() const {
  // This operator doesn't simulate any variables -- it only produces diagnostics.
  return oops::ObsVariables();
}

void ObsBackgroundErrorVertInterp::print(std::ostream & os) const {
  os << "ObsBackgroundErrorVertInterp: config = " << parameters_ << std::endl;
}

}  // namespace ufo
