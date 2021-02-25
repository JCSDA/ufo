/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/backgrounderrorvertinterp/ObsBackgroundErrorVertInterp.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/backgrounderrorvertinterp/ObsBackgroundErrorVertInterp.interface.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

static ObsOperatorMaker<ObsBackgroundErrorVertInterp> maker("BackgroundErrorVertInterp");

ObsBackgroundErrorVertInterp::ObsBackgroundErrorVertInterp(const ioda::ObsSpace & odb,
                                                           const eckit::Configuration & config)
  : ObsOperatorBase(odb, config),
    odb_(odb)
{
  oops::Log::trace() << "ObsBackgroundErrorVertInterp constructor entered" << std::endl;

  parameters_.validateAndDeserialize(config);

  requiredVars_.push_back(interpolationLevels());
  // simulateObs() may be asked to interpolate the background errors of any simulated variables.
  // We need to assume the worst, i.e. that we'll need to interpolate all of them.
  const oops::Variables &obsvars = odb.obsvariables();
  for (size_t ivar = 0; ivar < obsvars.size(); ++ivar)
    requiredVars_.push_back(obsvars[ivar] + "_background_error");

  oops::Log::trace() << "ObsBackgroundErrorVertInterp created" << std::endl;
}

ObsBackgroundErrorVertInterp::~ObsBackgroundErrorVertInterp() {
  oops::Log::trace() << "ObsBackgroundErrorVertInterp destructed" << std::endl;
}

void ObsBackgroundErrorVertInterp::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                               ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "ObsBackgroundErrorVertInterp: simulateObs entered" << std::endl;

  const std::string &verticalCoordinate = parameters_.verticalCoordinate;
  const std::string levels = interpolationLevels();

  oops::Variables variables;
  if (parameters_.variables.value() != boost::none)
    for (const Variable &variable : *parameters_.variables.value())
      variables += variable.toOopsVariables();
  else
    variables = odb_.obsvariables();

  ufo_backgrounderrorvertinterp_fillobsdiags_f90(verticalCoordinate.size(),
                                                 verticalCoordinate.c_str(),
                                                 levels.size(),
                                                 levels.c_str(),
                                                 geovals.toFortran(), odb_, hofx.nlocs(),
                                                 variables,
                                                 ydiags.toFortran());

  oops::Log::trace() << "ObsBackgroundErrorVertInterp: simulateObs exit" <<  std::endl;
}

const oops::Variables & ObsBackgroundErrorVertInterp::requiredVars() const {
  return requiredVars_;
}

oops::Variables ObsBackgroundErrorVertInterp::simulatedVars() const {
  // This operator doesn't simulate any variables -- it only produces diagnostics.
  return oops::Variables();
}

void ObsBackgroundErrorVertInterp::print(std::ostream & os) const {
  os << "ObsBackgroundErrorVertInterp: config = " << parameters_ << std::endl;
}

std::string ObsBackgroundErrorVertInterp::interpolationLevels() const {
  return parameters_.interpolationLevels.value().value_or(
        parameters_.verticalCoordinate.value());
}

}  // namespace ufo
