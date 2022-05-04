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
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
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
  const std::string &sensor = options_.sensor.value();

  // Currently the code is designed only for SSMIS brightness temperatures from
  // channels 12 through 18, but a different sensor could use a different list of
  // channels and frequencies requiring a different block of input checks.
  if (sensor == "SSMIS") {
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

  } else if (sensor == "AMSUA" || sensor == "ATMS") {
    ASSERT(options_.ch238d.value() != boost::none && options_.ch314d.value() != boost::none);
    channels_ = {options_.ch238d.value().get(), options_.ch314d.value().get()};
    ASSERT(options_.ch238d.value().get() != 0 && options_.ch314d.value().get() != 0 &&
           channels_.size() == 2);

  } else {
    std::string errString = "Currently only SSMIS, AMSUA, ATMS sensor are supported.";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }

  // required variables
  if (sensor == "AMSUA" || sensor == "ATMS") {
    geovars_ += oops::Variables({"water_area_fraction",
                                 "average_surface_temperature_within_field_of_view"});
    hdiags_ += oops::Variables({"brightness_temperature"}, vars.channels());
  }
}

// -----------------------------------------------------------------------------

void CloudLiquidWater::compute(const ioda::ObsSpace & odb,
                             const GeoVaLs & geovals,
                             const ObsDiagnostics & ydiags,
                             const ObsBias & biascoeffs,
                             ioda::ObsVector & out) const {
  // Get required parameters
  const std::string &vargrp = options_.varGroup.value();
  const std::string &sensor = options_.sensor.value();
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();

  const float fmiss = util::missingValue(fmiss);
  const double dmiss = util::missingValue(dmiss);

  std::vector<float> bt19h, bt19v, bt22v, bt37h, bt37v, bt91v, bt91h;
  if (sensor == "SSMIS") {
    bt19h.resize(nlocs);
    bt19v.resize(nlocs);
    bt22v.resize(nlocs);
    bt37h.resize(nlocs);
    bt37v.resize(nlocs);
    bt91h.resize(nlocs);
    bt91v.resize(nlocs);

    // From the obs database, grab all the brightness temperatures of channels.
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[0]), bt19h);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[1]), bt19v);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[2]), bt22v);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[3]), bt37h);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[4]), bt37v);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[5]), bt91v);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[6]), bt91h);
  }

  std::vector<float> bt238o, bt314o, bt238f, bt314f, bt238fBC, bt314fBC;
  if (sensor == "AMSUA" || sensor == "ATMS") {
    bt238o.resize(nlocs);
    bt314o.resize(nlocs);
    bt238f.resize(nlocs);
    bt314f.resize(nlocs);
    bt238fBC.resize(nlocs);
    bt314fBC.resize(nlocs);

    // From the obs database, grab all the brightness temperatures of channels.
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[0]), bt238o);
    odb.get_db(vargrp, "brightness_temperature_" + std::to_string(channels_[1]), bt314o);

    std::string hdiags;
    hdiags = "brightness_temperature_" + std::to_string(channels_[0]);
    ydiags.get(bt238f, hdiags);
    hdiags = "brightness_temperature_" + std::to_string(channels_[1]);
    ydiags.get(bt314f, hdiags);

    std::vector<float> scanangle(nlocs);
    odb.get_db("MetaData", "sensor_view_angle", scanangle);

    const Predictors & predictors = biascoeffs.predictors();
    const std::size_t npreds = predictors.size();
    double beta1, beta2;
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      beta1 = biascoeffs(0, channels_[0]-1);
      beta2 = biascoeffs(0, channels_[1]-1);
      bt238fBC[jloc] = bt238f[jloc] + beta1;
      bt314fBC[jloc] = bt314f[jloc] + beta2;
      for (std::size_t jord = 0; jord < 4; ++jord) {
        beta1 = biascoeffs(npreds-jord-1, channels_[0]-1);
        beta2 = biascoeffs(npreds-jord-1, channels_[1]-1);
        bt238fBC[jloc] += beta1 * pow(scanangle[jloc] * Constants::deg2rad, jord+1);
        bt314fBC[jloc] += beta2 * pow(scanangle[jloc] * Constants::deg2rad, jord+1);
      }
    }
  }

  std::vector<float> szas(nlocs);
  odb.get_db("MetaData", "sensor_zenith_angle", szas);

  std::vector<float> water_frac(nlocs, 0.0);
  std::vector<float> tsavg(nlocs, 0.0);
  geovals.get(water_frac, "water_area_fraction");
  if (sensor == "AMSUA" || sensor == "ATMS") {
    geovals.get(tsavg, "average_surface_temperature_within_field_of_view");
  }

  // Compute cloud liquid water amount
  std::vector<float> clw(nlocs, 0.0);
  if (sensor == "SSMIS") {
    CLWRetMW_SSMIS::cloudLiquidWater(bt19h, bt19v, bt22v, bt37h, bt37v, bt91v, bt91h,
                                     water_frac, clw);
  }
  if (sensor == "AMSUA" || sensor == "ATMS") {
    CloudLiquidWater::clwDerivative_amsua(tsavg, water_frac, bt238o, bt314o,
                                          bt238fBC, bt314fBC, clw);
  }

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
void CloudLiquidWater::clwDerivative_amsua(const std::vector<float> & tsavg,
                                           const std::vector<float> & water_frac,
                                           const std::vector<float> & bt238o,
                                           const std::vector<float> & bt314o,
                                           const std::vector<float> & bt238fBC,
                                           const std::vector<float> & bt314fBC,
                                           std::vector<float> & out) {
  ///
  /// \brief derivative of cloud liquid water from AMSU-A 23.8 GHz and 31.4 GHz channels.
  ///
  /// Reference: Grody et al. (2001)
  /// Determination of precipitable water and cloud liquid water over oceans from
  /// the NOAA 15 advanced microwave sounding unit
  ///
  const float t0c = Constants::t0c;
  const float d1 = 0.754, d2 = -2.265;
  const float c1 = 8.240, c2 = 2.622, c3 = 1.846;
  const float tbmax = 550.0, r284 = 284.0, r285 = 285.0;
  const float fmiss = util::missingValue(fmiss);

  for (size_t iloc = 0; iloc < water_frac.size(); ++iloc) {
    if (water_frac[iloc] >= 0.99) {
      if (tsavg[iloc] > t0c - 1.0 && bt238o[iloc] < tbmax && bt314o[iloc] < tbmax
                                  && bt238o[iloc] > 0.0 && bt314o[iloc] > 0.0
                                  && bt238fBC[iloc] <= r284 && bt314fBC[iloc] <= r284) {
        out[iloc] = d1 * (bt238fBC[iloc]-bt238o[iloc])/(r285 - bt238fBC[iloc])
                    + d2 * (bt314fBC[iloc]-bt314o[iloc])/(r285 - bt314fBC[iloc]);
      } else {
        out[iloc] = fmiss;
      }
    }
  }
}
// -----------------------------------------------------------------------------

}  // namespace ufo
