/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/backgrounderroridentity/ObsBackgroundErrorIdentity.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/backgrounderroridentity/ObsBackgroundErrorIdentity.interface.h"

namespace ufo {

static ObsOperatorMaker<ObsBackgroundErrorIdentity> maker("BackgroundErrorIdentity");

ObsBackgroundErrorIdentity::ObsBackgroundErrorIdentity(const ioda::ObsSpace & odb,
                                                       const Parameters_ & parameters)
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
    odb_(odb), parameters_(parameters)
{
  oops::Log::trace() << "ObsBackgroundErrorIdentity constructor entered" << std::endl;

  // simulateObs() may be asked to interpolate the background errors of any simulated variables.
  // We need to assume the worst, i.e. that we'll need to interpolate all of them.
  const oops::ObsVariables &obsvars = odb.assimvariables();
  for (size_t ivar = 0; ivar < obsvars.size(); ++ivar)
    requiredVars_.push_back(nameMap_.convertName(obsvars[ivar]).name() + "_background_error");

  oops::Log::trace() << "ObsBackgroundErrorIdentity created" << std::endl;
}

ObsBackgroundErrorIdentity::~ObsBackgroundErrorIdentity() {
  oops::Log::trace() << "ObsBackgroundErrorIdentity destructed" << std::endl;
}

void ObsBackgroundErrorIdentity::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                             ObsDiagnostics & ydiags,
                                             const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsBackgroundErrorIdentity: simulateObs entered" << std::endl;

  oops::Variables variables;
  if (parameters_.variables.value() != boost::none)
    for (const Variable &variable : *parameters_.variables.value())
      variables += nameMap_.convertName(variable.toOopsObsVariables());
  else
    variables = nameMap_.convertName(odb_.assimvariables());

  ufo_backgrounderroridentity_fillobsdiags_f90(geovals.toFortran(), hofx.nlocs(), variables,
                                               ydiags.toFortran());

  oops::Log::trace() << "ObsBackgroundErrorIdentity: simulateObs exit" <<  std::endl;
}

const oops::Variables & ObsBackgroundErrorIdentity::requiredVars() const {
  return requiredVars_;
}

oops::ObsVariables ObsBackgroundErrorIdentity::simulatedVars() const {
  // This operator doesn't simulate any variables -- it only produces diagnostics.
  return oops::ObsVariables();
}

void ObsBackgroundErrorIdentity::print(std::ostream & os) const {
  os << "ObsBackgroundErrorIdentity: config = " << parameters_ << std::endl;
}

}  // namespace ufo
