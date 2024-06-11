/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/filters/obsfunctions/CLWRetMW_SSMIS.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/predictors/CloudLiquidWater.h"
#include "ufo/utils/Constants.h"
namespace ufo {

static PredictorMaker<CloudLiquidWater>
       makerFuncCloudLiquidWater_("cloudWaterContent");

// -----------------------------------------------------------------------------

CloudLiquidWater::CloudLiquidWater(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1)) {
  if (parameters.order.value() != boost::none) {
    // override the predictor name by "cloudWaterContent_order_2"
    name() = name() + "_order_" + std::to_string(order_);
  }

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

  } else if (sensor == "GMI_GPM" || sensor == "gmi_gpm") {
    ASSERT((options_.ch37v.value() != boost::none && options_.ch37h.value() != boost::none) ||
           (options_.clwret_ch37v.value() != boost::none &&
            options_.clwret_ch37h.value() != boost::none));
    if (options_.ch37v.value() != boost::none && options_.ch37h.value() != boost::none) {
      channels_ = {options_.ch37v.value().get(), options_.ch37h.value().get()};
    } else {
      channels_ = {options_.clwret_ch37v.value().get(), options_.clwret_ch37h.value().get()};
    }
    if (options_.order.value() != boost::none) {
      order_ = options_.order.value().get();
    }
    hdiags_ += oops::ObsVariables({"brightness_temperature"}, vars.channels());
    hdiags_ += oops::ObsVariables({"brightness_temperature_assuming_clear_sky"}, vars.channels());
    hdiags_ += oops::ObsVariables({"transmittances_of_atmosphere_layer"}, vars.channels());
    geovars_ += oops::Variables({"water_area_fraction",
                                 "air_temperature",
                                 "air_pressure",
                                 "average_surface_temperature_within_field_of_view"});
//  Read tlap rate from the "tlapse_file". Copied this part from "LapseRate.cc"
    const std::string & tlapse_file = options_.tlapse.value().get();
    std::ifstream infile(tlapse_file);
    std::string nusis;   //  sensor/instrument/satellite
    int nuchan;  //  channel number
    float tlapse;

    if (infile.is_open()) {
      while (!infile.eof()) {
        infile >> nusis;
        infile >> nuchan;
        infile >> tlapse;
        tlapmean_[nuchan] = tlapse;
      }
      infile.close();
    }
//  End of reading tlapmean.
  } else {
    std::string errString = "Currently only SSMIS, AMSUA, ATMS sensor are supported.";
    oops::Log::error() << errString;
    throw eckit::BadValue(errString);
  }

  // required variables
  if (sensor == "AMSUA" || sensor == "ATMS") {
    geovars_ += oops::Variables({oops::Variable{"water_area_fraction"},
                              oops::Variable{"average_surface_temperature_within_field_of_view"}});
    hdiags_ += oops::ObsVariables({"brightness_temperature"}, vars.channels());
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

  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();

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
    odb.get_db(vargrp, "brightnessTemperature", bt19h, {channels_[0]});
    odb.get_db(vargrp, "brightnessTemperature", bt19v, {channels_[1]});
    odb.get_db(vargrp, "brightnessTemperature", bt22v, {channels_[2]});
    odb.get_db(vargrp, "brightnessTemperature", bt37h, {channels_[3]});
    odb.get_db(vargrp, "brightnessTemperature", bt37v, {channels_[4]});
    odb.get_db(vargrp, "brightnessTemperature", bt91v, {channels_[5]});
    odb.get_db(vargrp, "brightnessTemperature", bt91h, {channels_[6]});
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
    odb.get_db(vargrp, "brightnessTemperature", bt238o, {channels_[0]});
    odb.get_db(vargrp, "brightnessTemperature", bt314o, {channels_[1]});

    std::string hdiags;
    hdiags = "brightness_temperature_" + std::to_string(channels_[0]);
    ydiags.get(bt238f, hdiags);
    hdiags = "brightness_temperature_" + std::to_string(channels_[1]);
    ydiags.get(bt314f, hdiags);

    std::vector<float> scanangle(nlocs);
    odb.get_db("MetaData", "sensorViewAngle", scanangle);

    ASSERT(biascoeffs.nrecs() == 1);

    const Predictors & predictors = biascoeffs.predictors();
    const std::size_t npreds = predictors.size();
    double beta1, beta2;
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      beta1 = biascoeffs(0, channels_[0]-1, 0);
      beta2 = biascoeffs(0, channels_[1]-1, 0);
      bt238fBC[jloc] = bt238f[jloc] + beta1;
      bt314fBC[jloc] = bt314f[jloc] + beta2;
      for (std::size_t jord = 0; jord < 4; ++jord) {
        beta1 = biascoeffs(0, channels_[0]-1, npreds-jord-1);
        beta2 = biascoeffs(0, channels_[1]-1, npreds-jord-1);
        bt238fBC[jloc] += beta1 * pow(scanangle[jloc] * Constants::deg2rad, jord+1);
        bt314fBC[jloc] += beta2 * pow(scanangle[jloc] * Constants::deg2rad, jord+1);
      }
    }
  }

  std::vector<float> szas(nlocs);
  odb.get_db("MetaData", "sensorZenithAngle", szas);

  std::vector<float> water_frac(nlocs, 0.0);
  std::vector<float> tsavg(nlocs, 0.0);
  geovals.get(water_frac, oops::Variable{"water_area_fraction"});
  if (sensor == "AMSUA" || sensor == "ATMS") {
    geovals.get(tsavg, oops::Variable{"average_surface_temperature_within_field_of_view"});
  }

  // Compute cloud liquid water amount
  std::vector<float> clw(nlocs, 0.0);
  std::vector<float> clw_gmi_ch1_4(nlocs, 0.0);
  if (sensor == "SSMIS") {
    CLWRetMW_SSMIS::cloudLiquidWater(bt19h, bt19v, bt22v, bt37h, bt37v, bt91v, bt91h,
                                     water_frac, clw);
  }
  if (sensor == "GMI_GPM" || sensor == "gmi_gpm") {
    // Note: brightness_temperature_ and brightness_temperature_assuming_clear_sky_ at 37V GHz
    // and 37H GHz channels have to be corrected with the "constant", "lapse_rate_order_2",
    // "lapse_rate", and "scan_angle_order_4", "scan_angle_order_3", "scan_angle_order_2",
    // and "scan_angle" bias correction terms before being used to retrieve cloud indices by
    // the function "CLWRetMW::CIret_37v37h_diff".

    // Indices of data at channeles 37v and 37h in the array "channels"
    const int jch37v = 0;
    const int jch37h = 1;

    bt37h.resize(nlocs);
    bt37v.resize(nlocs);
    odb.get_db("ObsValue", "brightnessTemperature", bt37v, {channels_[jch37v]});
    odb.get_db("ObsValue", "brightnessTemperature", bt37h, {channels_[jch37h]});

    std::vector<float> bt_hofx_37vo(nlocs), bt_hofx_37ho(nlocs);
    ydiags.get(bt_hofx_37vo, "brightness_temperature_" + std::to_string(channels_[jch37v]));
    ydiags.get(bt_hofx_37ho, "brightness_temperature_" + std::to_string(channels_[jch37h]));

    std::vector<float> bt_clr_37vo(nlocs), bt_clr_37ho(nlocs);
    ydiags.get(bt_clr_37vo, "brightness_temperature_assuming_clear_sky_" +
               std::to_string(channels_[jch37v]));
    ydiags.get(bt_clr_37ho, "brightness_temperature_assuming_clear_sky_" +
               std::to_string(channels_[jch37h]));

    // retrieve the average surface temperature
    std::vector<float> tsavg5(nlocs, 0.0);
    geovals.get(tsavg5, oops::Variable{"average_surface_temperature_within_field_of_view"});

    std::vector<int> scanpos(nlocs, 0);
    odb.get_db("MetaData", "sensorScanPosition", scanpos);

    // common vectors storage, ptau5
    std::vector <float> pred(nlocs, 0.0);
    // Retrieve the transmittances_of_atmosphere_layer from Hdiag
    std::vector<std::vector<std::vector<float>>> ptau5;
    std::vector<std::vector<float>> tmpvar;

    std::string hdiags;
    for (std::size_t jvar = 0; jvar < channels_.size(); ++jvar) {
      hdiags = "transmittances_of_atmosphere_layer_" + std::to_string(channels_[jvar]);
      tmpvar.clear();
      for (std::size_t js = 0; js < ydiags.nlevs(hdiags); ++js) {
        ydiags.get(pred, hdiags, js);
        tmpvar.push_back(pred);
      }
      ptau5.push_back(tmpvar);
    }

    // Retrieve the temperature,tvp
    std::vector<std::vector<float>> tvp;
    std::size_t nlevs = geovals.nlevs(oops::Variable{"air_temperature"});
    for (std::size_t js = 0; js < nlevs; ++js) {
      geovals.getAtLevel(pred, oops::Variable{"air_temperature"}, js);
      tvp.push_back(pred);
    }

    // sort out the tlapmean based on vars, tlap
    std::vector<float> tlap;
    for (std::size_t jvar = 0; jvar < channels_.size(); ++jvar) {
      auto it = tlapmean_.find(channels_[jvar]);
      if (it != tlapmean_.end()) {
        tlap.push_back(it->second);
      } else {
        oops::Log::error() << "Could not locate tlapemean for channel: " <<
                              channels_[jvar] << std::endl;
        ABORT("Could not locate tlapemean value");
      }
    }

    clw_bias_correction_gmi(biascoeffs,
                          bt_hofx_37vo, bt_hofx_37ho, bt_clr_37vo, bt_clr_37ho, bt37v, bt37h,
                           water_frac, tsavg5, scanpos, ptau5, tvp, tlap, channels_, nlevs,
                           clw, clw_gmi_ch1_4);
  }
  if (sensor == "AMSUA" || sensor == "ATMS") {
    CloudLiquidWater::clwDerivative_amsua(tsavg, water_frac, bt238o, bt314o,
                                          bt238fBC, bt314fBC, clw);
  }
  if (sensor != "GMI_GPM" && sensor != "gmi_gpm") {
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
  } else {
    for (std::size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        if (vars_.channels()[jvar] >= 4) {
          out[iloc*nvars + jvar] = static_cast<double>(pow(clw[iloc], order_));
        } else {
          out[iloc*nvars + jvar] = static_cast<double>(pow(clw_gmi_ch1_4[iloc], order_));
        }
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
  const float fmiss = util::missingValue<float>();

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

void CloudLiquidWater::clw_bias_correction_gmi(const ObsBias & biascoeffs,
             const std::vector<float> & bt_hofx_37vo,
             const std::vector<float> & bt_hofx_37ho,
             const std::vector<float> & bt_clr_37vo,
             const std::vector<float> & bt_clr_37ho,
             const std::vector<float> & bt37v,
             const std::vector<float> & bt37h,
             const std::vector<float> & water_frac,
             const std::vector<float> & tsavg5,
             const std::vector<int> & scanpos,
             const std::vector<std::vector<std::vector<float>>> & ptau5,
             const std::vector<std::vector<float>> & tvp,
             const std::vector<float> & tlap,
             const std::vector<int>  & channels_,
             const int & nlevs,
             std::vector<float> & clw,
             std::vector<float> & clw_gmi_ch1_4) {
    // Indices of data at channeles 37v and 37h in the array "channels"
    const int jch37v = 0;
    const int jch37h = 1;
    //  Beginning constant and scan_angle terms
    const std::size_t nlocs = bt_hofx_37vo.size();
    std::vector<float> bias_37v(nlocs), bias_37h(nlocs);
    const Predictors & predictors = biascoeffs.predictors();
    const std::size_t npreds = predictors.size();
    double beta1, beta2;

    std::vector<std::string> predictors_part = {"constant", "lapseRate_order_2", "lapseRate",
      "sensorScanAngle_order_4", "sensorScanAngle_order_3", "sensorScanAngle_order_2",
      "sensorScanAngle"};
    size_t id_preds = predictors_part.size();
    std::vector<int> id_pred(id_preds, -1);
    for (std::size_t kp = 0; kp < id_preds; ++kp) {
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        const std::string varname = predictors[jp]->name();
        if (varname == predictors_part[kp]) {
          id_pred[kp] = jp;
          break;
        }
      }
      ASSERT(id_pred[kp] >= 0);
    }

    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      // constant term
      beta1 = biascoeffs(0, channels_[jch37v]-1, id_pred[0]);
      beta2 = biascoeffs(0, channels_[jch37h]-1, id_pred[0]);
      bias_37v[jloc] = beta1;
      bias_37h[jloc] = beta2;
      // scan_position term
      for (std::size_t jord = 0; jord < 4; ++jord) {
        beta1 = biascoeffs(0, channels_[jch37v]-1, id_pred[id_preds-jord-1]);
        beta2 = biascoeffs(0, channels_[jch37h]-1, id_pred[id_preds-jord-1]);
        bias_37v[jloc] += beta1 * pow(scanpos[jloc] * Constants::deg2rad, jord+1);
        bias_37h[jloc] += beta2 * pow(scanpos[jloc] * Constants::deg2rad, jord+1);
      }
    }

//  Beginning of tlap BC for GMI data. Modified from "LapseRate.cc".
    float tlapchn;
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
//    For channel_37v
      int jvar = jch37v;
      tlapchn = (ptau5[jvar][nlevs-2][jloc]-ptau5[jvar][nlevs-1][jloc])*
                (tsavg5[jloc]-tvp[nlevs-2][jloc]);
      for (std::size_t k = 1; k < nlevs-1; ++k) {
        tlapchn = tlapchn+(ptau5[jvar][nlevs-k-2][jloc]-ptau5[jvar][nlevs-k-1][jloc])*
                  (tvp[nlevs-k][jloc]-tvp[nlevs-k-2][jloc]);
      }
      bias_37v[jloc] += (biascoeffs(0, channels_[jvar]-1, id_pred[1]) *
                        pow((tlapchn - tlap[jvar]), 2));
      bias_37v[jloc] += (biascoeffs(0, channels_[jvar]-1, id_pred[2])*(tlapchn - tlap[jvar]));
//    For channel_37h
      jvar = jch37h;
      tlapchn = (ptau5[jvar][nlevs-2][jloc]-ptau5[jvar][nlevs-1][jloc])*
                (tsavg5[jloc]-tvp[nlevs-2][jloc]);
      for (std::size_t k = 1; k < nlevs-1; ++k) {
        tlapchn = tlapchn+(ptau5[jvar][nlevs-k-2][jloc]-ptau5[jvar][nlevs-k-1][jloc])*
                  (tvp[nlevs-k][jloc]-tvp[nlevs-k-2][jloc]);
      }
      bias_37h[jloc] += (biascoeffs(0, channels_[jvar]-1, id_pred[1]) *
                        pow((tlapchn - tlap[jvar]), 2));
      bias_37h[jloc] += (biascoeffs(0, channels_[jvar]-1, id_pred[2])*(tlapchn - tlap[jvar]));
    }
//  End of tlap BC for GMI data

//  Sum up contant, tlap, and scan position BC for BC
    std::vector<float> bt_hofx_37v(nlocs), bt_hofx_37h(nlocs);
    std::vector<float> bt_clr_37v(nlocs), bt_clr_37h(nlocs);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      bt_hofx_37v[jloc] = bt_hofx_37vo[jloc] + bias_37v[jloc];
      bt_hofx_37h[jloc] = bt_hofx_37ho[jloc] + bias_37h[jloc];
      bt_clr_37v[jloc]  = bt_clr_37vo[jloc] + bias_37v[jloc];
      bt_clr_37h[jloc]  = bt_clr_37ho[jloc] + bias_37h[jloc];
    }

    // Retrieve cloud liquid water indices
    std::vector<float> clw_obs(nlocs), clw_hofx(nlocs);
    CLWRetMW::CIret_37v37h_diff(bt_clr_37v, bt_clr_37h, water_frac, bt37v, bt37h, clw_obs);
    CLWRetMW::CIret_37v37h_diff(bt_clr_37v, bt_clr_37h, water_frac,
                                bt_hofx_37v, bt_hofx_37h, clw_hofx);

    float clw_ret;
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      int cld_rbc_idx = 1;
//    Copied from GSI subroutine: "radiance_mod.f90, radiance_ex_biascor_gmi()"
      if (clw_obs[iloc] > 0.05 || clw_hofx[iloc] > 0.05 ||
          abs(clw_obs[iloc] - clw_hofx[iloc]) >= 0.001) {
        cld_rbc_idx = 0;
      }
//    For GMI channels 4-13
      if (clw_obs[iloc] > 0.05 && clw_hofx[iloc] > 0.05 &&
          clw_obs[iloc] != bad_clwret_value_  && clw_hofx[iloc] != bad_clwret_value_ &&
        cld_rbc_idx == 0) {
        clw_ret = 0.5 * (clw_obs[iloc] + clw_hofx[iloc]);
        clw[iloc] = clw_ret;
      }
//    For GMI channels 1-3)
      if (clw_obs[iloc] > 0.2 && clw_hofx[iloc] > 0.2 &&
          clw_obs[iloc] != bad_clwret_value_  && clw_hofx[iloc] != bad_clwret_value_ &&
        cld_rbc_idx == 0) {
        clw_ret = 0.5 * (clw_obs[iloc] + clw_hofx[iloc]);
        clw_gmi_ch1_4[iloc] = clw_ret;
      }
    }
  }
// -----------------------------------------------------------------------------

}  // namespace ufo
