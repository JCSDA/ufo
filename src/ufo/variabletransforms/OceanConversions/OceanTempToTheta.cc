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

  // First copy across the ObsError information to derived variable if present,
  // by doing everything before adding the derived obs we will set
  // QCflag::missing based on the presence or absence of theta rather than the
  // auxillary variables.
  {
    std::vector<float> thetaError;
    getObservation("ObsError", temperaturevariable_, thetaError);
    const size_t iv = obserr_.varnames().find(thetavariable_);
    if (!thetaError.empty()) {
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (!apply[jobs])
          continue;
        obserr_[iv][jobs] = thetaError[jobs];
      }
      putObservation(thetavariable_, thetaError, "DerivedObsError");
    }
  }

  // Copy across the GrossErrorProbability information to derived variable if present
  {
    std::vector<float> thetapge;
    getObservation("GrossErrorProbability", temperaturevariable_,
                   thetapge);
    if (!thetapge.empty())
      putObservation(thetavariable_, thetapge, "GrossErrorProbability");
  }

  // Copy across QC flag information
  {
    std::vector<int> thetaflags;
    getObservation("QCflagsData", temperaturevariable_,
                   thetaflags);
    if (thetaflags.empty()) {
      thetaflags.assign(nlocs, 0);
    }
    putObservation(thetavariable_, thetaflags, "QCflagsData");
  }

  // compute theta as function of temperature, pressure and salinity
  {
    std::vector<float> pressure;
    std::vector<float> temp;
    std::vector<float> sal;
    std::vector<float> theta(nlocs, missingValueFloat);
    getObservation(salinitygroup_,
                   salinityvariable_,
                   sal, true);
    getObservation(temperaturegroup_,
                   temperaturevariable_,
                   temp, true);
    getObservation(pressuregroup_,
                   pressurevariable_,
                   pressure, true);
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
    putObservation(thetavariable_, theta, "DerivedObsValue");
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
