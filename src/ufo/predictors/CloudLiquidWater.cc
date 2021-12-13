/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/CLWRetMW_SSMIS.h"
#include "ufo/GeoVaLs.h"
#include "ufo/predictors/CloudLiquidWater.h"
#include "ufo/utils/Constants.h"
namespace ufo {

static PredictorMaker<CloudLiquidWater>
       makerFuncCloudLiquidWater_("cloud_liquid_water");

// -----------------------------------------------------------------------------

CloudLiquidWater::CloudLiquidWater(const Parameters_ & parameters, const oops::Variables & vars)
  : PredictorBase(parameters, vars) {
  // Initialize options
  options_ = parameters;
  const std::string &satellite = options_.satellite.value();

  // Currently the code is designed only for SSMIS brightness temperatures from
  // channels 12 through 18, but a different satellite could use a different list of
  // channels and frequencies requiring a different block of input checks.
  if (satellite == "SSMIS") {
    ASSERT(options_.ch19h.value() != boost::none && options_.ch19v.value() != boost::none &&
           options_.ch22v.value() != boost::none && options_.ch37h.value() != boost::none &&
           options_.ch37v.value() != boost::none && options_.ch91v.value() != boost::none &&
           options_.ch91h.value() != boost::none);

    channels_ = {options_.ch19h.value().get(), options_.ch19v.value().get(),
                 options_.ch22v.value().get(), options_.ch37h.value().get(),
                 options_.ch37v.value().get(), options_.ch91v.value().get(),
                 options_.ch91h.value().get()};

    ASSERT(options_.ch19h.value().get() != 0 && options_.ch19v.value().get() != 0 &&
           options_.ch22v.value().get() != 0 && options_.ch37h.value().get() != 0 &&
           options_.ch37v.value().get() != 0 && options_.ch91v.value().get() != 0 &&
           options_.ch91h.value().get() != 0 && channels_.size() == 7);

  } else {
    std::string errString = "Currently only SSMIS satellite is supported.";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }
}

// -----------------------------------------------------------------------------

void CloudLiquidWater::compute(const ioda::ObsSpace & odb,
                             const GeoVaLs & geovals,
                             const ObsDiagnostics &,
                             ioda::ObsVector & out) const {
  // Get required parameters
  const std::string &vargrp = options_.varGroup.value();
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();

  const float fmiss = util::missingValue(fmiss);
  const double dmiss = util::missingValue(dmiss);

  // From the obs database, grab all the brightness temperatures of channels.
  std::vector<float> bt19h(nlocs), bt19v(nlocs), bt22v(nlocs), bt37h(nlocs),
                     bt37v(nlocs), bt91v(nlocs), bt91h(nlocs);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[0]), bt19h);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[1]), bt19v);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[2]), bt22v);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[3]), bt37h);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[4]), bt37v);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[5]), bt91v);
  odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[6]), bt91h);

  std::vector<float> szas(nlocs);
  odb.get_db("MetaData", "sensor_zenith_angle", szas);

  std::vector<float> water_frac(nlocs, 0.0);
  geovals.get(water_frac, "water_area_fraction");

  // Compute cloud liquid water amount
  std::vector<float> clw(nlocs);
  CLWRetMW_SSMIS::cloudLiquidWater(bt19h, bt19v, bt22v, bt37h, bt37v, bt91v, bt91h,
                                   water_frac, clw);

  for (std::size_t iloc = 0; iloc < nlocs; ++iloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      float cossza = cos(Constants::deg2rad * szas[iloc]);
      if (clw[iloc] == fmiss) {
        out[iloc*nvars + jvar] = dmiss;
      } else {
        out[iloc*nvars + jvar] = static_cast<double>(clw[iloc]*cossza*cossza);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
