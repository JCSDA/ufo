/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CloudDetectMinResidualAVHRR.h"

#include <cmath>

#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<CloudDetectMinResidualAVHRR>
       makerCloudDetectMinResidualAVHRR_("CloudDetectMinResidualAVHRR");

// -----------------------------------------------------------------------------

CloudDetectMinResidualAVHRR::CloudDetectMinResidualAVHRR(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Include required variables from ObsDiag
  invars_ += Variable("ObsDiag/brightness_temperature_jacobian_surface_temperature", channels_);
  invars_ += Variable("ObsDiag/brightness_temperature_jacobian_air_temperature", channels_);
  invars_ += Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_);
  invars_ += Variable("ObsDiag/pressure_level_at_peak_of_weightingfunction", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);
  invars_ += Variable(errgrp+"/brightnessTemperature", channels_);
  invars_ += Variable(hofxgrp+"/brightnessTemperature", channels_);
  invars_ += Variable("ObsValue/brightnessTemperature", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
  invars_ += Variable("GeoVaLs/land_area_fraction");
  invars_ += Variable("GeoVaLs/ice_area_fraction");
  invars_ += Variable("GeoVaLs/surface_snow_area_fraction");
  invars_ += Variable("GeoVaLs/average_surface_temperature_within_field_of_view");
  invars_ += Variable("GeoVaLs/air_pressure");
  invars_ += Variable("GeoVaLs/air_temperature");
  invars_ += Variable("GeoVaLs/tropopause_pressure");
}

// -----------------------------------------------------------------------------

CloudDetectMinResidualAVHRR::~CloudDetectMinResidualAVHRR() {}

// -----------------------------------------------------------------------------

void CloudDetectMinResidualAVHRR::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "CloudDetectMinResidualAVHRR::compute start" << std::endl;
  // Get channel use flags from options
  std::vector<int> use_flag = options_.useflagChannel.value();
  std::vector<int> use_flag_clddet = options_.useflagCloudDetect.value();

  // Get tuning parameters for surface sensitivity over sea/land/oce/snow/mixed from options
  std::vector<float> dtempf_in = options_.obserrScaleFactorTsfc.value();

  // Get dimensions
  size_t nlocs = in.nlocs();
  if (nlocs == 0) {
    return;
  }
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("GeoVaLs/air_pressure"));

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Get variables from ObsDiag
  // Load surface temperature jacobian
  std::vector<std::vector<float>> dbtdts(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsDiag/brightness_temperature_jacobian_surface_temperature",
                     channels_)[ichan], dbtdts[ichan]);
  }

  // Get temperature jacobian
  std::vector<std::vector<std::vector<float>>>
       dbtdt(nchans, std::vector<std::vector<float>>(nlevs, std::vector<float>(nlocs)));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      const size_t level = nlevs - ilev - 1;
      in.get(Variable("ObsDiag/brightness_temperature_jacobian_air_temperature",
                       channels_)[ichan], level, dbtdt[ichan][ilev]);
    }
  }

  // Get layer-to-space transmittance
  std::vector<std::vector<std::vector<float>>>
       tao(nchans, std::vector<std::vector<float>>(nlevs, std::vector<float>(nlocs)));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      const size_t level = nlevs - ilev - 1;
      in.get(Variable("ObsDiag/transmittances_of_atmosphere_layer",
             channels_)[ichan], level,  tao[ichan][ilev]);
    }
  }

  // Get pressure level at the peak of the weighting function
  std::vector<float> values(nlocs, 0.0);
  std::vector<std::vector<float>> wfunc_pmaxlev(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsDiag/pressure_level_at_peak_of_weightingfunction",
                     channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      wfunc_pmaxlev[ichan][iloc] = nlevs - values[iloc] + 1;
    }
  }

  // Get variables from ObsSpace
  // Get effective observation error and convert it to inverse of the error variance
  const float missing = util::missingValue<float>();
  std::vector<int> qcflag(nlocs, 0);
  std::vector<std::vector<float>> varinv_use(nchans, std::vector<float>(nlocs, 0.0));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], values);
    in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflag);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (flaggrp == "PreQC") values[iloc] == missing ? qcflag[iloc] = 100 : qcflag[iloc] = 0;
      (qcflag[iloc] == 0) ? (values[iloc] = 1.0 / pow(values[iloc], 2)) : (values[iloc] = 0.0);
      if (use_flag_clddet[ichan] > 0 && use_flag_clddet[ichan]%2 == 1)
          varinv_use[ichan][iloc] = values[iloc];
    }
  }

  // Get bias corrected innovation (tbobs - hofx - bias)
  std::vector<std::vector<float>> innovation(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsValue/brightnessTemperature", channels_)[ichan], innovation[ichan]);
    in.get(Variable(hofxgrp+"/brightnessTemperature", channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innovation[ichan][iloc] = innovation[ichan][iloc] - values[iloc];
    }
  }

  // Get original observation error.  If not explicitly passed through the YAML, check the ObsSpace.
  // This channel-dependent error is assumed to be constant across all obs locations.
  std::vector<float> obserr(nchans, 0.0f);
  if (options_.obserrOriginal.value() != boost::none) {
    obserr = options_.obserrOriginal.value().get();
  } else {
  // Get original observation error (uninflated) from ObsSpace
      std::vector<std::vector<float>> obserr2(nchans, std::vector<float>(nlocs));
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      in.get(Variable("ObsError/brightnessTemperature", channels_)[ichan], obserr2[ichan]);
      obserr[ichan] = obserr2[ichan][0];
    }
  }

  // Get tropopause pressure [Pa]
  std::vector<float> tropprs(nlocs);
  in.get(Variable("GeoVaLs/tropopause_pressure"), tropprs);

  // Get average surface temperature within FOV
  std::vector<float> tsavg(nlocs);
  in.get(Variable("GeoVaLs/average_surface_temperature_within_field_of_view"), tsavg);

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  std::vector<float> land_frac(nlocs);
  std::vector<float> ice_frac(nlocs);
  std::vector<float> snow_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);
  in.get(Variable("GeoVaLs/land_area_fraction"), land_frac);
  in.get(Variable("GeoVaLs/ice_area_fraction"), ice_frac);
  in.get(Variable("GeoVaLs/surface_snow_area_fraction"), snow_frac);

  bool sea = true;
  bool land = false;
  bool ice = false;
  bool snow = false;
  bool mixed = false;
  float dtempf;

  // Get air pressure [Pa]
  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    in.get(Variable("GeoVaLs/air_pressure"), level, prsl[ilev]);
  }

  // Get air temperature
  std::vector<std::vector<float>> tair(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    in.get(Variable("GeoVaLs/air_temperature"), level, tair[ilev]);
  }

  // Minimum Residual Method (MRM) for Cloud Detection:
  // Determine model level index of the cloud top (lcloud)
  // Find pressure of the cloud top (cldprs)
  // Estimate cloud fraction (cldfrac)
  // output: out = 0 clear channel
  //         out = 1 cloudy channel
  //         out = 2 clear channel but too sensitive to surface condition

  // Loop through locations
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    float sum3 = 0.0;
    float tmp = 0.0;
    float cloudp = 0.0;
    sea = water_frac[iloc] >= 0.99;
    land = land_frac[iloc] >= 0.99;
    ice = ice_frac[iloc] >= 0.99;
    snow = snow_frac[iloc] >= 0.99;
    mixed = (!sea && !land && !ice && !snow);
    if (sea) {
      dtempf = dtempf_in[0];
    } else if (land) {
      dtempf = dtempf_in[1];
    } else if (ice) {
      dtempf = dtempf_in[2];
    } else if (snow) {
      dtempf = dtempf_in[3];
    } else {
      dtempf = dtempf_in[4];
    }
    std::vector<std::vector<float>> dbt(nchans, std::vector<float>(nlevs));
    for (size_t ichan=0; ichan < nchans; ++ichan) {
      if (varinv_use[ichan][iloc] > 0) {
        sum3 = sum3 + innovation[ichan][iloc] * innovation[ichan][iloc] * varinv_use[ichan][iloc];
      }
    }
    sum3 = 0.75 * sum3;

    // Set initial cloud condition
    int lcloud = 0;
    float cldfrac = 0.0;
    float cldprs = prsl[0][iloc] * 0.01;     // convert from [Pa] to [hPa]

    // Loop through vertical layer from surface to model top
    for (size_t k = 0 ; k < nlevs ; ++k) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) dbt[ichan][k] = 0.0;
      // Perform cloud detection within troposphere
      if (prsl[k][iloc] * 0.01 > tropprs[iloc] * 0.01) {
        float sum = 0.0, sum2 = 0.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          if (varinv_use[ichan][iloc] > 0.0) {
            dbt[ichan][k] = (tair[k][iloc] - tsavg[iloc]) * dbtdts[ichan][iloc];
            for (size_t kk = 0; kk < k; ++kk) {
              dbt[ichan][k] = dbt[ichan][k] + (tair[k][iloc] - tair[kk][iloc]) *
                              dbtdt[ichan][kk][iloc];
            }
            sum = sum + innovation[ichan][iloc] * dbt[ichan][k] * varinv_use[ichan][iloc];
            sum2 = sum2 +  dbt[ichan][k] * dbt[ichan][k] * varinv_use[ichan][iloc];
          }
        }
        if (fabs(sum2) < FLT_MIN) sum2 = copysign(1.0e-12, sum2);
        cloudp = std::min(std::max((sum/sum2), 0.f), 1.f);
        sum = 0.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          if (varinv_use[ichan][iloc] > 0.0) {
          tmp = innovation[ichan][iloc] - cloudp * dbt[ichan][k];
          sum = sum + tmp * tmp * varinv_use[ichan][iloc];
          }
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
    // Cloud check
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      size_t ilev;
      out[ichan][iloc] = 0;
      for (ilev = 0; ilev < lcloud; ++ilev) {
        if (fabs(cldfrac * dbt[ichan][ilev]) > obserr[ichan]) {
          out[ichan][iloc]= 1;
          varinv_use[ichan][iloc]= 0.0;
          break;
        }
      }
    }
    // If no clouds is detected, do sensivity to surface temperature check
    // Initialize at each location
    float sumx = 0.0, sumx2 = 0.0;
    float dts = 0.0;
    float delta = 0.0;
    const float dts_threshold = 3.0;
    for (size_t ichan=0; ichan < nchans; ++ichan) {
      sumx = sumx + innovation[ichan][iloc] * dbtdts[ichan][iloc] * varinv_use[ichan][iloc];
      sumx2 = sumx2 + dbtdts[ichan][iloc] * dbtdts[ichan][iloc] * varinv_use[ichan][iloc];
    }
    if (fabs(sumx2) < FLT_MIN) sumx2 = copysign(1.0e-12, sumx2);
    dts = std::fabs(sumx / sumx2);
    float dts_save = dts;
    if (std::abs(dts) > 1.0) {
      if (sea == false) {
        dts = std::min(dtempf, dts);
      } else {
        dts = std::min(dts_threshold, dts);
      }
      for (size_t ichan=0; ichan < nchans; ++ichan) {
        delta = obserr[ichan];
        if (std::abs(dts * dbtdts[ichan][iloc]) > delta) out[ichan][iloc] = 2;
      }
    }
  // end of location loop
  }
  oops::Log::trace() << "CloudDetectMinResidualAVHRR::compute done" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & CloudDetectMinResidualAVHRR::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
