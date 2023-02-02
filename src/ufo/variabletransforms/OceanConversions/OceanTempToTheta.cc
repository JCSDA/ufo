/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/OceanConversions/OceanTempToTheta.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<OceanTempToTheta>
    makerOceanTempToTheta_("OceanTempToTheta");

OceanTempToTheta::OceanTempToTheta(
        const Parameters_ &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      salinityvariable_(options.SalinityVariable),
      salinitygroup_(options.SalinityGroup),
      temperaturevariable_(options.TemperatureVariable),
      temperaturegroup_(options.TemperatureGroup),
      pressurevariable_(options.PressureVariable),
      pressuregroup_(options.PressureGroup),
      thetavariable_(options.ThetaVariable)
{}

// -----------------------------------------------------------------------------

void OceanTempToTheta::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() <<
    "Retrieve Ocean Potential Temperature from Salinity, Temperature, Pressure" <<
    std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Get all required data
  std::vector<float> pressure;
  std::vector<float> temp;
  std::vector<float> sal;
  std::vector<float> theta(nlocs, missingValueFloat);
  std::vector<float> thetaError(nlocs, missingValueFloat);
  std::vector<float> thetapge;
  std::vector<int> thetaflags;
  getObservation(salinitygroup_,
                 salinityvariable_,
                 sal, true);
  getObservation(temperaturegroup_,
                 temperaturevariable_,
                 temp, true);
  getObservation(pressuregroup_,
                 pressurevariable_,
                 pressure, true);
  data_.get(Variable(std::string("ObsErrorData/") + temperaturevariable_), thetaError);
  getObservation("GrossErrorProbability", temperaturevariable_,
                 thetapge);
  getObservation("QCFlags", temperaturevariable_,
                 thetaflags);

  if (thetapge.empty()) {
    thetapge.assign(nlocs, missingValueFloat);
  }
  if (thetaflags.empty()) {
    thetaflags.assign(nlocs, 0);
  }

  // compute theta as function of temperature, pressure and salinity
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (sal[loc] != missingValueFloat &&
        temp[loc] != missingValueFloat &&
        pressure[loc] != missingValueFloat) {
      theta[loc] = gsw_pt_from_t_f90(sal[loc],
                                     temp[loc],
                                     pressure[loc]);
    }
  }
  obsdb_.put_db("DerivedObsValue", thetavariable_, theta);
  const size_t iv = obserr_.varnames().find(thetavariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = thetaError[jobs];
  }

  // copy ObsError, PGEFinal and QCflags to new potential temperature
  obsdb_.put_db("GrossErrorProbability", thetavariable_, thetapge);
  obsdb_.put_db("QCFlags", thetavariable_, thetaflags);
  obsdb_.put_db("DerivedObsError", thetavariable_, thetaError);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
