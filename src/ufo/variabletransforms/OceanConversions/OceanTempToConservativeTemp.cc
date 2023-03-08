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
{
    if (!obserr_.varnames().has(temperaturevariable_) ||
        !obserr_.varnames().has(conservativetempvariable_) ||
        !flags_.varnames().has(temperaturevariable_) ||
        !flags_.varnames().has(conservativetempvariable_)) {
        throw eckit::BadValue("Both `" + conservativetempvariable_ + "` and `"
            + temperaturevariable_ + "` must be observed or derived variables for" +
            " the `OceanTempToConservativeTemp` variable transform.", Here());
    }
}

// -----------------------------------------------------------------------------

void OceanTempToConservativeTemp::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() <<
    "Retrieve Ocean Conservative Temperature from Absolute Salinity, In Situ Temp., Pressure" <<
    std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Copy across the ObsErrorData and QCflagsData information to derived
  // variable if present, by doing everything before adding the derived obs we
  // will set QCflag::missing based on the presence or absence of the
  // observation rather than the auxillary variables.
  {
    const size_t iErrTemp = obserr_.varnames().find(temperaturevariable_);
    const size_t iErrConTemp = obserr_.varnames().find(conservativetempvariable_);
    const size_t iFlagTemp = flags_.varnames().find(temperaturevariable_);
    const size_t iFlagConTemp = flags_.varnames().find(conservativetempvariable_);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!apply[jobs])
        continue;
      obserr_[iErrConTemp][jobs] = obserr_[iErrTemp][jobs];
      flags_[iFlagConTemp][jobs] = flags_[iFlagTemp][jobs];
    }
  }

  // Copy across the GrossErrorProbability information to derived variable if present
  {
    std::vector<float> conservativetemppge;
    getObservation("GrossErrorProbability", temperaturevariable_,
                   conservativetemppge);
    if (!conservativetemppge.empty())
      putObservation(conservativetempvariable_, conservativetemppge, "GrossErrorProbability");
  }

  // compute conservative temp as function of temperature, pressure and salinity
  {
    std::vector<float> pressure;
    std::vector<float> temp;
    std::vector<float> sal;
    std::vector<float> conservativetemp(nlocs, missingValueFloat);
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
        conservativetemp[loc] = gsw_ct_from_t_f90(sal[loc],
                                                  temp[loc],
                                                  pressure[loc]);
      }
    }
    putObservation(conservativetempvariable_, conservativetemp, "DerivedObsValue");
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
