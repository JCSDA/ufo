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
{}

// -----------------------------------------------------------------------------

void OceanPracticalSalinityToAbsoluteSalinity::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() <<
    "Retrieve Ocean Absolute Salinity from Practical Salinity, Pressure, Latitude and Longitude" <<
    std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Get all required data
  std::vector<float> pressure;
  std::vector<float> psal;
  std::vector<float> lats;
  std::vector<float> lons;
  std::vector<float> asal(nlocs, missingValueFloat);
  std::vector<float> asalError(nlocs, missingValueFloat);
  std::vector<float> asalpge;
  std::vector<int> asalflags;
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
  data_.get(Variable(std::string("ObsErrorData/") + practicalsalinityvariable_), asalError);
  getObservation("GrossErrorProbability", practicalsalinityvariable_,
                 asalpge);
  getObservation("QCFlags", practicalsalinityvariable_,
                 asalflags);

  if (asalpge.empty()) {
    asalpge.assign(nlocs, missingValueFloat);
  }
  if (asalflags.empty()) {
    asalflags.assign(nlocs, 0);
  }

  // compute absolute salinity as a function of practical salinity, pressure, longitude and latitude
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
  obsdb_.put_db("DerivedObsValue", absolutesalinityvariable_, asal);
  const size_t iv = obserr_.varnames().find(absolutesalinityvariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = asalError[jobs];
  }

  // copy ObsError, PGEFinal and QCflags to new conservative temperature
  obsdb_.put_db("GrossErrorProbability", absolutesalinityvariable_, asalpge);
  obsdb_.put_db("QCFlags", absolutesalinityvariable_, asalflags);
  obsdb_.put_db("DerivedObsError", absolutesalinityvariable_, asalError);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
