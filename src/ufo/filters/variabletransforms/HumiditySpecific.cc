/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/variabletransforms/HumiditySpecific.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/variabletransforms/HumidityRelative.h"
#include "ufo/utils/Constants.h"

namespace ufo {

/* -------------------------------------------------------------------------------------*/

HumiditySpecific::HumiditySpecific(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::trace() << "HumiditySpecific constructor starting" << std::endl;
  options_.deserialize(config);
  std::string methodStr = options_.method_for_RH.value();
  oops::Log::debug() << "    method for RH is: " << methodStr << std::endl;
}

/* -------------------------------------------------------------------------------------*/

HumiditySpecific::~HumiditySpecific() {}

/* -------------------------------------------------------------------------------------*/

void HumiditySpecific::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "HumiditySpecific applyFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();
  std::string errString;

  // Get the method for calculating RH from the options.
  std::string methodStr = options_.method_for_RH.value();
  auto method = formulasRH::resolveMethods(methodStr);

  std::vector<float> rh(nlocs), t(nlocs), p(nlocs);
  std::vector<float> specific_humidity(nlocs);
  // output is assumed missing unless a valid value calculated below
  specific_humidity.assign(nlocs, missing);

  if (obsdb_.has("ObsValue", "relative_humidity")) {
    obsdb_.get_db("ObsValue", "relative_humidity", rh);
  } else {
    errString = "Unable to continue because relative_humidity@ObsValue is missing";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }
  if (obsdb_.has("ObsValue", "air_temperature")) {
    obsdb_.get_db("ObsValue", "air_temperature", t);
  } else {
    errString = "Unable to continue because air_temperature@ObsValue is missing";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }
  if (obsdb_.has("MetaData", "air_pressure")) {
    obsdb_.get_db("MetaData", "air_pressure", p);
  } else if (obsdb_.has("ObsValue", "surface_pressure")) {
    obsdb_.get_db("ObsValue", "surface_pressure", p);
  } else {
    errString = "Unable to continue, requires either air_pressure@MetaData or surface_pressure";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }

  float esat, qvs, qv;

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // Check for missing values.
      if (rh[jobs] != missing && t[jobs] != missing && p[jobs] != missing) {
        // Cannot realistically continue with pathological inputs.
        if (t[jobs] > 180 && t[jobs] < 350 && p[jobs] > 0.0001 && p[jobs] < 110000) {
          // Calculate saturation vapor pressure from temperature according to requested method
          // Double-check result is always lower than 15% of incoming pressure.
          esat = std::min(p[jobs]*0.15f, formulasRH::SatVaporPres_fromTemp(t[jobs], method));
          // Convert sat. vapor pressure to sat water vapor mixing ratio
          qvs = 0.622 * esat/(p[jobs]-esat);
          // Using RH, calculate water vapor mixing ratio
          qv = std::max(1.0e-12f, rh[jobs]*qvs);
          // Final specific humidity is kept to a sensible lowest limit
          specific_humidity[jobs] = std::max(1.0e-12f, qv/(1.0f+qv));
        }
      }
    }
  }
  // put new variable at existing locations using DerivedValue group
  if (options_.save) obsdb_.put_db("DerivedValue", "specific_humidity", specific_humidity);
}

/* -------------------------------------------------------------------------------------*/

void HumiditySpecific::print(std::ostream & os) const {
  os << "HumiditySpecific";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
