/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/variabletransforms/HumidityRelative.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace formulasRH {

/* -------------------------------------------------------------------------------------*/

Method resolveMethods(const std::string & input) {
  if (input == "UKMO") return UKMO;
  if (input == "NOAA") return NOAA;
  if (input == "Walko") return Walko;
  if (input == "Murphy") return Murphy;
  return EMPTY;
}

/* -------------------------------------------------------------------------------------*/

float SatVaporPres_fromTemp(float temp_K, Method method) {
  float e_sub_s;  // Saturation vapor pressure (over water) output.
  const float t0 = static_cast<float>(ufo::Constants::t0c);

  switch (method) {
    case formulasRH::Method::UKMO: {
    /* Source: Eqn 7, Sonntag, D., Advancements in the field of hygrometry,
     *                Meteorol. Zeitschrift, N. F., 3, 51-66, 1994.
     * Most radiosonde manufacturers use Wexler, or Hyland and Wexler
     * or Sonntag formulations, which are all very similar (Holger Vomel, pers. comm., 2011)
    */
      e_sub_s = std::exp(-6096.9385f / temp_K + 21.2409642f - 2.711193E-2f * temp_K +
                      1.673952E-5f * temp_K * temp_K + 2.433502f * std::log(temp_K));
      break;
    }

    case formulasRH::Method::Walko: {
      // Polynomial fit of Goff-Gratch (1946) formulation. (Walko, 1991)
      float x = std::max(-80.0f, temp_K-t0);
      const std::vector<float> c{610.5851f, 44.40316f, 1.430341f, 0.2641412e-1f,
                   0.2995057e-3f, 0.2031998e-5f, 0.6936113e-8f, 0.2564861e-11f, -0.3704404e-13f};
      e_sub_s = c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*(c[4]+x*(c[5]+x*(c[6]+x*(c[7]+x*c[8])))))));
      break;
    }

    case formulasRH::Method::Murphy: {
      // ALTERNATIVE (costs more CPU, more accurate than Walko, 1991)
      // Source: Murphy and Koop, Review of the vapour pressure of ice and
      //       supercooled water for atmospheric applications, Q. J. R.
      //       Meteorol. Soc (2005), 131, pp. 1539-1565.

      e_sub_s = std::exp(54.842763f - 6763.22f / temp_K - 4.210f * std::log(temp_K)
                     + 0.000367f * temp_K + std::tanh(0.0415f * (temp_K - 218.8f))
                     * (53.878f - 1331.22f / temp_K - 9.44523f * std::log(temp_K)
                     + 0.014025f * temp_K));
      break;
    }

    case formulasRH::Method::EMPTY: {
      std::string errString = "Aborting, no method matches enum formulasRH::Method";
      oops::Log::error() << errString;
      throw eckit::BadValue(errString);
      break;
    }

    case formulasRH::Method::NOAA:
    default: {
      // Classical formula from Rogers and Yau (1989; Eq2.17)
      e_sub_s = 1000.*0.6112*std::exp(17.67f*(temp_K-t0)/(temp_K-29.65f));
      break;
    }
  }
  return e_sub_s;
}

}  // namespace formulasRH

namespace ufo {

/* -------------------------------------------------------------------------------------*/

HumidityRelative::HumidityRelative(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::trace() << "HumidityRelative constructor starting" << std::endl;
  options_.deserialize(config);
  std::string methodStr = options_.method_for_RH.value();
  oops::Log::debug() << "    method for RH is: " << methodStr << std::endl;
}

/* -------------------------------------------------------------------------------------*/

HumidityRelative::~HumidityRelative() {}

/* -------------------------------------------------------------------------------------*/

void HumidityRelative::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "HumidityRelative applyFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();
  std::string errString;

  // Get the method for calculating RH from the options.
  std::string methodStr = options_.method_for_RH.value();
  auto method = formulasRH::resolveMethods(methodStr);

  std::vector<float> sh(nlocs), t(nlocs), p(nlocs);
  std::vector<float> relative_humidity(nlocs);
  // output is assumed missing unless a valid value calculated below
  relative_humidity.assign(nlocs, missing);

  if (obsdb_.has("ObsValue", "specific_humidity")) {
    obsdb_.get_db("ObsValue", "specific_humidity", sh);
  } else {
    errString = "Unable to continue because specific_humidity@ObsValue is missing";
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
      if (sh[jobs] != missing && t[jobs] != missing && p[jobs] != missing) {
        // Cannot realistically continue with pathological inputs.
        sh[jobs] = std::max(1.0e-9f, sh[jobs]);
        if (t[jobs] > 180 && t[jobs] < 350 && p[jobs] > 0.0001 && p[jobs] < 110000) {
          // Calculate saturation vapor pressure from temperature according to requested method
          // Double-check result is always lower than 15% of incoming pressure.
          esat = std::min(p[jobs]*0.15f, SatVaporPres_fromTemp(t[jobs], method));
          // Convert sat. vapor pressure to sat water vapor mixing ratio
          qvs = 0.622 * esat/(p[jobs]-esat);
          // Convert specific humidity to water vapor mixing ratio
          qv = std::max(1.0e-12f, sh[jobs]/(1.0f-sh[jobs]));
          // Final RH (which can be greater than 100%) is q/qsat, but set sensible lowest limit
          relative_humidity[jobs] = std::max(1.0e-6f, qv/qvs);
        }
      }
    }
  }
  // put new variable at existing locations using DerivedValue group
  if (options_.save) obsdb_.put_db("DerivedValue", "relative_humidity", relative_humidity);
}

/* -------------------------------------------------------------------------------------*/

void HumidityRelative::print(std::ostream & os) const {
  os << "HumidityRelative";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
