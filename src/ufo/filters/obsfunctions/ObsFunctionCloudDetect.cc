/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionCloudDetect.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionCloudDetect> makerObsFuncCloudDetect_("CloudDetect");

// -----------------------------------------------------------------------------

ObsFunctionCloudDetect::ObsFunctionCloudDetect(const eckit::LocalConfiguration conf)
  : invars_(), errgrp_("ObsErrorData"), hofxgrp_("HofX"), biasgrp_("ObsBias"),
    channels_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("channels") && conf_.has("use_flag") && conf_.has("use_flag_clddet") &&
         conf_.has("obserr_dtempf"));

  // Check if using user-defined obserr (should be available from obs file) for testing
  if (conf_.has("obserr_test")) errgrp_ = conf_.getString("obserr_test");
  if (conf_.has("hofx_test")) hofxgrp_ = conf_.getString("hofx_test");
  if (conf_.has("bias_test")) biasgrp_ = conf_.getString("bias_test");

  // Get channels from options
  const std::string chlist = conf.getString("channels");
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include required variables from ObsDiag
  invars_ += Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag", channels_);
  invars_ += Variable("brightness_temperature_jacobian_air_temperature@ObsDiag", channels_);
  invars_ += Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_);
  invars_ += Variable("pressure_level_at_peak_of_weightingfunction@ObsDiag", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@"+errgrp_, channels_);
  invars_ += Variable("brightness_temperature@"+biasgrp_, channels_);
  invars_ += Variable("brightness_temperature@"+hofxgrp_, channels_);
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@ObsError", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("land_area_fraction@GeoVaLs");
  invars_ += Variable("ice_area_fraction@GeoVaLs");
  invars_ += Variable("surface_snow_area_fraction@GeoVaLs");
  invars_ += Variable("average_surface_temperature_within_field_of_view@GeoVaLs");
  invars_ += Variable("air_pressure@GeoVaLs");
  invars_ += Variable("air_temperature@GeoVaLs");
  invars_ += Variable("tropopause_pressure@GeoVaLs");
}

// -----------------------------------------------------------------------------

ObsFunctionCloudDetect::~ObsFunctionCloudDetect() {}

// -----------------------------------------------------------------------------

void ObsFunctionCloudDetect::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get channel use flags from options
  std::vector<int> use_flag = conf_.getIntVector("use_flag");
  std::vector<int> use_flag_clddet = conf_.getIntVector("use_flag_clddet");

  // Get tuning parameters for surface sensitivity over sea/land/oce/snow/mixed from options
  std::vector<float> dtempf_in = conf_.getFloatVector("obserr_dtempf");

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("air_pressure@GeoVaLs"));

  // Get variables from ObsDiag
  // Load surface temperature jacobian
  std::vector<std::vector<float>> dbtdts(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag",
                     channels_)[ichan], dbtdts[ichan]);
  }

  // Get temperature jacobian
  std::vector<std::vector<std::vector<float>>>
       dbtdt(nchans, std::vector<std::vector<float>>(nlevs, std::vector<float>(nlocs)));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      int level = nlevs - ilev;
      in.get(Variable("brightness_temperature_jacobian_air_temperature@ObsDiag",
                       channels_)[ichan], level, dbtdt[ichan][ilev]);
    }
  }

  // Get layer-to-space transmittance
  std::vector<std::vector<std::vector<float>>>
       tao(nchans, std::vector<std::vector<float>>(nlevs, std::vector<float>(nlocs)));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      int level = nlevs - ilev;
      in.get(Variable("transmittances_of_atmosphere_layer@ObsDiag",
             channels_)[ichan], level,  tao[ichan][ilev]);
    }
  }

  // Get pressure level at the peak of the weighting function
  std::vector<float> values(nlocs, 0.0);
  std::vector<std::vector<float>> wfunc_pmaxlev(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("pressure_level_at_peak_of_weightingfunction@ObsDiag",
                     channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      wfunc_pmaxlev[ichan][iloc] = nlevs - values[iloc] + 1;
    }
  }

  // Get variables from ObsSpace
  // Get effective observation error and convert it to inverse of the error variance
  std::vector<std::vector<float>> varinv_use(nchans, std::vector<float>(nlocs, 0.0));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+errgrp_, channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      values[iloc] = 1.0 / pow(values[iloc], 2);
      if (use_flag_clddet[ichan] > 0) varinv_use[ichan][iloc] = values[iloc];
    }
  }

  // Get bias corrected innovation (tbobs - hofx - bias)
  std::vector<std::vector<float>> innovation(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ichan], innovation[ichan]);
    in.get(Variable("brightness_temperature@"+hofxgrp_, channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innovation[ichan][iloc] = innovation[ichan][iloc] - values[iloc];
    }
    in.get(Variable("brightness_temperature@"+biasgrp_, channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innovation[ichan][iloc] = innovation[ichan][iloc] - values[iloc];
    }
  }

  // Get original observation error (uninflated) from ObsSpaec
  std::vector<std::vector<float>> obserr(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsError", channels_)[ichan], obserr[ichan]);
  }

  // Get variables from GeoVaLS
  // Get tropopause pressure [Pa]
  std::vector<float> tropprs(nlocs);
  in.get(Variable("tropopause_pressure@GeoVaLs"), tropprs);

  // Get average surface temperature within FOV
  std::vector<float> tsavg(nlocs);
  in.get(Variable("average_surface_temperature_within_field_of_view@GeoVaLs"), tsavg);

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  std::vector<float> land_frac(nlocs);
  std::vector<float> ice_frac(nlocs);
  std::vector<float> snow_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);
  in.get(Variable("land_area_fraction@GeoVaLs"), land_frac);
  in.get(Variable("ice_area_fraction@GeoVaLs"), ice_frac);
  in.get(Variable("surface_snow_area_fraction@GeoVaLs"), snow_frac);

  // Determine dominant surface type in each FOV
  std::vector<bool> land(nlocs, false);
  std::vector<bool> sea(nlocs, false);
  std::vector<bool> ice(nlocs, false);
  std::vector<bool> snow(nlocs, false);
  std::vector<bool> mixed(nlocs, false);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    sea[iloc] = water_frac[iloc] >= 0.99;
    land[iloc] = land_frac[iloc] >= 0.99;
    ice[iloc] = ice_frac[iloc] >= 0.99;
    snow[iloc] = snow_frac[iloc] >= 0.99;
    mixed[iloc] = (!sea[iloc] && !land[iloc] && !ice[iloc] && !snow[iloc]);
  }

  // Setup weight given to each surface type
  std::vector<float> dtempf(nlocs);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (sea[iloc]) {
      dtempf[iloc] = dtempf_in[0];
    } else if (land[iloc]) {
      dtempf[iloc] = dtempf_in[1];
    } else if (ice[iloc]) {
      dtempf[iloc] = dtempf_in[2];
    } else if (snow[iloc]) {
      dtempf[iloc] = dtempf_in[3];
    } else {
      dtempf[iloc] = dtempf_in[4];
    }
  }

  // Get air pressure [Pa]
  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    size_t level = nlevs - ilev;
    in.get(Variable("air_pressure@GeoVaLs"), level, prsl[ilev]);
  }

  // Get air temperature
  std::vector<std::vector<float>> tair(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    size_t level = nlevs - ilev;
    in.get(Variable("air_temperature@GeoVaLs"), level, tair[ilev]);
  }

  // Minimum Residual Method (MRM) for Cloud Detection:
  // Determine model level index of the cloud top (lcloud)
  // Find pressure of the cloud top (cldprs)
  // Estimate cloud fraction (cldfrac)
  // output: out = 0 clear channel
  //         out = 1 cloudy channel
  //         out = 2 clear channel but too sensitive to surface condition

  // Set vectors to hold cloud infomation from cloud detection (cab be part of the output)
  // std::vector<float> cloud_prsl(nlocs);
  // std::vector<float> cloud_frac(nlocs);
  // std::vector<int> cloud_lev(nlocs);

  // Loop through locations
  const float btmax = 550.0, btmin = 50.0;
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    // Initialize at each location
    // cloud_lev[iloc] = 0;
    // cloud_prsl[iloc] = 0.0;
    // cloud_frac[iloc] = 0.0;
    float sum = 0.0, sum2 = 0.0, sum3 = 0.0;
    float tmp = 0.0;
    float cloudp = 0.0;
    std::vector<float> dbt(nchans);
    for (size_t ichan=0; ichan < nchans; ++ichan) {
      sum3 = sum3 + innovation[ichan][iloc] * innovation[ichan][iloc] * varinv_use[ichan][iloc];
    }
    sum3 = 0.75 * sum3;

    // Set initial cloud condition
    int lcloud = 0;
    float cldfrac = 0.0;
    float cldprs = prsl[0][iloc] * 0.01;     // convert from [Pa] to [hPa]

    // Loop through vertical layer from surface to model top
    for (size_t k = 0 ; k < nlevs ; ++k) {
      // Perform cloud detection within troposphere
      if (prsl[k][iloc] * 0.01 > tropprs[iloc] * 0.01) {
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          dbt[ichan] = (tair[k][iloc] - tsavg[iloc]) * dbtdts[ichan][iloc];
        }
        for (size_t kk = 0; kk < k; ++kk) {
          for (size_t ichan = 0; ichan < nchans; ++ichan) {
            dbt[ichan] = dbt[ichan] + (tair[k][iloc] - tair[kk][iloc]) * dbtdt[ichan][kk][iloc];
          }
        }
        sum = 0.0;
        sum2 = 0.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          if (varinv_use[ichan][iloc] > 0.0) {
            sum = sum + innovation[ichan][iloc] * dbt[ichan] * varinv_use[ichan][iloc];
            sum2 = sum2 +  dbt[ichan] * dbt[ichan] * varinv_use[ichan][iloc];
          }
        }
        cloudp = std::min(std::max(static_cast<double>(sum/sum2), 0.0), 1.0);
        sum = 0.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
         // if (varinv_use[ichan][iloc] > 0.0) {
          tmp = innovation[ichan][iloc] - cloudp * dbt[ichan];
          sum = sum + tmp * tmp * varinv_use[ichan][iloc];
         // }
        }
        if (sum < sum3) {
          sum3 = sum;
          lcloud = k + 1;   // array index + 1 -> model coordinate index
          cldfrac = cloudp;
          cldprs = prsl[k][iloc] * 0.01;
        }
      }
    // end of vertical loop
    }
    // If more than 2% of the transmittance comes from the cloud layer,
    // reject the channel (marked as cloudy channel)
    float tao_cld = -999.0;
    for (size_t ichan = 0; ichan < nchans; ++ichan) out[ichan][iloc] = 0;
    if (lcloud > 0) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        // Get cloud top transmittance
        tao_cld = tao[ichan][lcloud-1][iloc];
        // Passive channels
        if (use_flag[ichan] < 0 && lcloud  >= wfunc_pmaxlev[ichan][iloc]) out[ichan][iloc] = 1;
        // Active channels
        if (out[ichan][iloc] != 1 && tao_cld > 0.02) out[ichan][iloc] = 1;
      }
      // cloud infomation output at model level
      // cloud_lev[iloc] = lcloud;
      // cloud_prsl[iloc] = cldprs;
      // cloud_frac[iloc] = cldfrac;
    } else {
    // If no clouds is detected, do sensivity to surface temperature check
    // Initialize at each location
      float sum = 0.0, sum2 = 0.0;
      float dts = 0.0;
      float delta = 0.0;
      const float dts_threshold = 3.0;
      for (size_t ichan=0; ichan < nchans; ++ichan) {
        delta = 0.0;
        sum = sum + innovation[ichan][iloc] * dbtdts[ichan][iloc] * varinv_use[ichan][iloc];
        sum2 = sum2 + dbtdts[ichan][iloc] * dbtdts[ichan][iloc] * varinv_use[ichan][iloc];
      }
      dts = std::fabs(sum / sum2);
      if (std::abs(dts) > 1.0) {
        if (sea[iloc] == false) {
          dts = std::min(dtempf[iloc], dts);
        } else {
          dts = std::min(dts_threshold, dts);
        }
        for (size_t ichan=0; ichan < nchans; ++ichan) {
          delta = std::max(0.05 * obserr[ichan][iloc], 0.02);
          if (std::abs(dts * dbtdts[ichan][iloc]) > delta) out[ichan][iloc] = 2;
        }
      }
    }
  // end of location loop
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionCloudDetect::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
