/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/OceanConversions/OceanTempToConservativeTemp.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<OceanTempToConservativeTemp>
    makerOceanTempToConservativeTemp_("OceanTempToConservativeTemp");

OceanTempToConservativeTemp::OceanTempToConservativeTemp(
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
      conservativetempvariable_(options.ConservativeTempVariable)
{}

// -----------------------------------------------------------------------------

void OceanTempToConservativeTemp::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() <<
    "Retrieve Ocean Conservative Temperature from Absolute Salinity, In Situ Temp., Pressure" <<
    std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Get all required data
  std::vector<float> pressure;
  std::vector<float> temp;
  std::vector<float> sal;
  std::vector<float> conservativetemp(nlocs, missingValueFloat);
  std::vector<float> conservativetempError(nlocs, missingValueFloat);
  std::vector<float> conservativetemppge;
  std::vector<int> conservativetempflags;
  getObservation(salinitygroup_,
                 salinityvariable_,
                 sal, true);
  getObservation(temperaturegroup_,
                 temperaturevariable_,
                 temp, true);
  getObservation(pressuregroup_,
                 pressurevariable_,
                 pressure, true);
  data_.get(Variable(std::string("ObsErrorData/") + temperaturevariable_), conservativetempError);
  getObservation("GrossErrorProbability", temperaturevariable_,
                 conservativetemppge);
  getObservation("QCFlags", temperaturevariable_,
                 conservativetempflags);

  if (conservativetemppge.empty()) {
    conservativetemppge.assign(nlocs, missingValueFloat);
  }
  if (conservativetempflags.empty()) {
    conservativetempflags.assign(nlocs, 0);
  }

  // compute conservative temp as function of temperature, pressure and salinity
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (sal[loc] != missingValueFloat &&
        temp[loc] != missingValueFloat &&
        pressure[loc] != missingValueFloat) {
      conservativetemp[loc] = gsw_ct_from_t_f90(sal[loc],
                                                temp[loc],
                                                pressure[loc]);
    }
  }
  obsdb_.put_db("DerivedObsValue", conservativetempvariable_, conservativetemp);
  const size_t iv = obserr_.varnames().find(conservativetempvariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = conservativetempError[jobs];
  }

  // copy ObsError, PGEFinal and QCflags to new conservative temperature
  obsdb_.put_db("GrossErrorProbability", conservativetempvariable_, conservativetemppge);
  obsdb_.put_db("QCFlags", conservativetempvariable_, conservativetempflags);
  obsdb_.put_db("DerivedObsError", conservativetempvariable_, conservativetempError);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
