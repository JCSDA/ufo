/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/OceanConversions/OceanPracticalSalinityToAbsoluteSalinity.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<OceanPracticalSalinityToAbsoluteSalinity>
    makerOceanPracticalSalinityToAbsoluteSalinity_("OceanPracticalSalinityToAbsoluteSalinity");

OceanPracticalSalinityToAbsoluteSalinity::OceanPracticalSalinityToAbsoluteSalinity(
        const Parameters_ &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      absolutesalinityvariable_(options.AbsoluteSalinityVariable),
      practicalsalinityvariable_(options.PracticalSalinityVariable),
      practicalsalinitygroup_(options.PracticalSalinityGroup),
      pressurevariable_(options.PressureVariable),
      pressuregroup_(options.PressureGroup)
{
    if (!obserr_.varnames().has(practicalsalinityvariable_) ||
        !obserr_.varnames().has(absolutesalinityvariable_) ||
        !flags_.varnames().has(practicalsalinityvariable_) ||
        !flags_.varnames().has(absolutesalinityvariable_)) {
        throw eckit::BadValue("Both `" + absolutesalinityvariable_ + "` and `"
            + practicalsalinityvariable_ + "` must be observed or derived variables" +
            " for the `OceanPracticalSalinityToAbsoluteSalinity` variable transform.", Here());
    }
}

// -----------------------------------------------------------------------------

void OceanPracticalSalinityToAbsoluteSalinity::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() <<
    "Retrieve Ocean Absolute Salinity from Practical Salinity, Pressure, Latitude and Longitude" <<
    std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Copy across observation error and QC flag info
  {
    const size_t iErrASal = obserr_.varnames().find(absolutesalinityvariable_);
    const size_t iErrPSal = obserr_.varnames().find(practicalsalinityvariable_);
    const size_t iFlagASal = flags_.varnames().find(absolutesalinityvariable_);
    const size_t iFlagPSal = flags_.varnames().find(practicalsalinityvariable_);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!apply[jobs])
        continue;
      obserr_[iErrASal][jobs] = obserr_[iErrPSal][jobs];
      flags_[iFlagASal][jobs] = flags_[iFlagPSal][jobs];
    }
  }

  // Copy across GrossErrorProbability if present
  {
    std::vector<float> asalpge;
    getObservation("GrossErrorProbability", practicalsalinityvariable_,
                   asalpge);
    if (!asalpge.empty())
      putObservation(absolutesalinityvariable_, asalpge, "GrossErrorProbability");
  }

  // compute absolute salinity as a function of practical salinity, pressure, longitude and latitude
  {
    std::vector<float> pressure;
    std::vector<float> psal;
    std::vector<float> lats;
    std::vector<float> lons;
    std::vector<float> asal(nlocs, missingValueFloat);
    getObservation(practicalsalinitygroup_,
                   practicalsalinityvariable_,
                   psal, true);
    getObservation(pressuregroup_,
                   pressurevariable_,
                   pressure, true);
    getObservation("MetaData",
                   "latitude",
                   lats, true);
    getObservation("MetaData",
                   "longitude",
                   lons, true);
    for (size_t loc = 0; loc < nlocs; ++loc) {
      if (!apply[loc]) continue;
      if (psal[loc] != missingValueFloat &&
          pressure[loc] != missingValueFloat) {
        asal[loc] = gsw_sa_from_sp_f90(psal[loc],
                                       pressure[loc],
                                       lons[loc],
                                       lats[loc]);
      }
    }
    putObservation(absolutesalinityvariable_, asal, "DerivedObsValue");
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
