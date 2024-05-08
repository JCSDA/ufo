/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/HydrometeorCheckAMSUAclr.h"

#include <cmath>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

static ObsFunctionMaker<HydrometeorCheckAMSUAclr>
                       makerHydrometeorCheckAMSUAclr_("HydrometeorCheckAMSUAclr");

// -----------------------------------------------------------------------------

HydrometeorCheckAMSUAclr::HydrometeorCheckAMSUAclr(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &biaspredgrp = options_.testBiasPred.value().get();
  const std::string &biastermgrp = options_.testBiasTerm.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Include required variables from ObsDiag
  invars_ += Variable("ObsDiag/brightness_temperature_jacobian_surface_emissivity", channels_);
  invars_ += Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable("ObsValue/brightnessTemperature", channels_);
  invars_ += Variable(hofxgrp+"/brightnessTemperature", channels_);
  invars_ += Variable(biastermgrp+"/constant"+biastermgrp, channels_);
  invars_ += Variable(biastermgrp+"/sensorScanAngle_order_4", channels_);
  invars_ += Variable(biastermgrp+"/sensorScanAngle_order_3", channels_);
  invars_ += Variable(biastermgrp+"/sensorScanAngle_order_2", channels_);
  invars_ += Variable(biastermgrp+"/sensorScanAngle", channels_);
  invars_ += Variable(biaspredgrp+"/brightnessTemperature", channels_);
  invars_ += Variable("MetaData/sensorZenithAngle");

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
  invars_ += Variable("GeoVaLs/land_area_fraction");
}

// -----------------------------------------------------------------------------

HydrometeorCheckAMSUAclr::~HydrometeorCheckAMSUAclr() {}

// -----------------------------------------------------------------------------

void HydrometeorCheckAMSUAclr::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);


  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Set channel index
  int ich238 = 0, ich314 = 1, ich503 = 2, ich528 = 3, ich536 = 4;
  int ich544 = 5, ich549 = 6, ich890 = 14;
  if (inst == "atms") {
    ich238 = 0, ich314 = 1, ich503 = 2, ich528 = 4, ich536 = 5;
    ich544 = 6, ich549 = 7, ich890 = 15;
  }

  // Get test groups from options
  const std::string &biaspredgrp = options_.testBiasPred.value().get();
  const std::string &biastermgrp = options_.testBiasTerm.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  std::vector<float> land_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);
  in.get(Variable("GeoVaLs/land_area_fraction"), land_frac);

  // Get surface temperature jacobian
  std::vector<std::vector<float>> dbtde(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsDiag/brightness_temperature_jacobian_surface_emissivity", channels_)[ichan],
           dbtde[ichan]);
  }

  // Get satellite zenith angle
  std::vector<float> satzen(nlocs);
  in.get(Variable("MetaData/sensorZenithAngle"), satzen);

  // Get HofX for clear-sky simulation
  std::vector<float> hofxclr536(nlocs);
  in.get(Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_)[ich536],
         hofxclr536);

  // Get bias predictor
  std::vector<float> clw_pred(nlocs);
  in.get(Variable(biaspredgrp+"/brightnessTemperature", channels_)[ich238], clw_pred);

  // Get ObsBiasTerm: constant term for 23.8GHz channel
  std::vector<float> bias_const238(nlocs);
  in.get(Variable(biastermgrp+"/constant", channels_)[ich238], bias_const238);

  // Get ObsBiasTerm: scan angle terms for 23.8GHz channel
  size_t nangs = 4;
  std::vector<float> values(nlocs);
  std::vector<std::string> scanterms(nangs);
  std::vector<float> bias_scanang238(nlocs, 0.0);
  scanterms[0] = biastermgrp+"/sensorScanAngle_order_4";
  scanterms[1] = biastermgrp+"/sensorScanAngle_order_3";
  scanterms[2] = biastermgrp+"/sensorScanAngle_order_2";
  scanterms[3] = biastermgrp+"/sensorScanAngle";
  for (size_t iang = 0; iang < nangs; ++iang) {
    in.get(Variable(scanterms[iang], channels_)[ich238], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      bias_scanang238[iloc] = bias_scanang238[iloc] + values[iloc];
    }
  }

  // Calculate bias-corrected innovation: Observation - HofX (HofX includes bias correction)
  std::vector<std::vector<float>> btobs(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> innov(nchans, std::vector<float>(nlocs));
  std::vector<float> hofx(nlocs);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsValue/brightnessTemperature", channels_)[ichan], btobs[ichan]);
    in.get(Variable(hofxgrp+"/brightnessTemperature", channels_)[ichan], hofx);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innov[ichan][iloc] = btobs[ichan][iloc] - hofx[iloc];
    }
  }

  // Set parameters
  const float fmissing = util::missingValue<float>();
  float w1f6 = 1.0/10.0, w2f6 = 1.0/0.80;
  float w1f4 = 1.0/0.30, w2f4 = 1.0/1.80;

  std::vector<std::vector<int>> affected_channels(nchans, std::vector<int>(nlocs));
  // Loop over locations
  // Combined cloud-precipitation-surface checks
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // Initialization and check on failed CLW calculation
    for (size_t ich = 0; ich < nchans; ++ich) {
      affected_channels[ich][iloc] = 0;
      if (clw_pred[iloc] == fmissing && ((ich <= ich544) || (ich >= ich890))) {
        affected_channels[ich][iloc] = 1;
      }
    }

    // Window channel sanity check
    // If any of the window channels is bad, skip all window channels
    // List of surface sensitivity channels
    std::vector<float> OmFs{std::abs(innov[ich238][iloc]),
                      std::abs(innov[ich314][iloc]), std::abs(innov[ich503][iloc]),
                      std::abs(innov[ich528][iloc]), std::abs(innov[ich536][iloc]),
                      std::abs(innov[ich544][iloc]), std::abs(innov[ich890][iloc])};
    bool result = false;
    result = any_of(OmFs.begin(), OmFs.end(), [](float x){
             return (x > 200.0 || x == util::missingValue<float>());});
    if (result) {
      for (size_t ich = 0; ich < nchans; ++ich) {
        if (affected_channels[ich][iloc] == 0 && (ich <= ich544 || ich >= ich890)) {
              affected_channels[ich][iloc] = 1;
        }
      }
    }

    if (!result && clw_pred[iloc] != fmissing) {
      // Calculate cloud liquid water and scattering effect
      float clwx = 0.6;
      float dsval = 0.8;
      if (water_frac[iloc] >= 0.99) {
        float cossza = cos(Constants::deg2rad * satzen[iloc]);
        clwx = w1f4 * clw_pred[iloc] / cossza;
        float btobsbc238 = btobs[ich238][iloc] - bias_const238[iloc]
                                               - bias_scanang238[iloc];
        dsval = ((2.410 - 0.0098 * btobsbc238) * innov[ich238][iloc] +
                0.454 * innov[ich314][iloc] - innov[ich890][iloc]) * w1f6;
        dsval = std::max(static_cast<float>(0.0), dsval);
      }
      float factch4 = pow(clwx, 2) + pow(innov[ich528][iloc] * w2f4, 2);
      float factch6 = pow(dsval, 2) + pow(innov[ich544][iloc] * w2f6, 2);

      // Hydrometeor check
      // Precipitation check (factch6)
      bool sea = water_frac[iloc] >= 0.99;
      bool land = land_frac[iloc] >= 0.99;
      bool latms_surfaceqc = (inst == "atms") && !(sea || land);

      if (factch6 >= 1.0 || latms_surfaceqc) {
        for (size_t ich = ich238; ich <= ich544; ++ich) {
          affected_channels[ich][iloc] = 1;
        }
        for (size_t ich = ich890; ich < nchans; ++ich) {
          affected_channels[ich][iloc] = 1;
        }
      // Thick cloud check (factch4)
      } else if (factch4 > 0.5) {
        for (size_t ich = ich238; ich <= ich536; ++ich) {
          affected_channels[ich][iloc] = 1;
        }
        for (size_t ich = ich890; ich < nchans; ++ich) {
          affected_channels[ich][iloc] = 1;
        }
      // Sensitivity of BT to the surface emissivity check
      } else {
        float thrd238, thrd314, thrd503, thrd890;
        if (water_frac[iloc] >= 0.99) {
          thrd238 = 0.025, thrd314 = 0.015, thrd503 = 0.030, thrd890 = 0.030;
        } else {
          thrd238 = 0.020, thrd314 = 0.015, thrd503 = 0.035, thrd890 = 0.015;
        }
        float de238 = 0.0, de314 = 0.0, de503 = 0.0, de890 = 0.0;
        float dbtde238 = dbtde[ich238][iloc];
        float dbtde314 = dbtde[ich314][iloc];
        float dbtde503 = dbtde[ich503][iloc];
        float dbtde890 = dbtde[ich890][iloc];
        if (dbtde238 != 0.0) de238 = std::abs(innov[ich238][iloc]) / dbtde238;
        if (dbtde314 != 0.0) de314 = std::abs(innov[ich314][iloc]) / dbtde314;
        if (dbtde503 != 0.0) de503 = std::abs(innov[ich503][iloc]) / dbtde503;
        if (dbtde890 != 0.0) de890 = std::abs(innov[ich890][iloc]) / dbtde890;
        bool qcemiss = false;
        qcemiss = de238 > thrd238 || de314 > thrd314 || de503 > thrd503 || de890 > thrd890;
        if (qcemiss) {
          for (size_t ich = ich238; ich <= ich536; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          for (size_t ich = ich890; ich < nchans; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
        }
      // end Hydrometeor check
      }
    // if over Window channel sanity check
    }
  // loop over locations
  }
  // Output
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ichan][iloc] = affected_channels[ichan][iloc];
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & HydrometeorCheckAMSUAclr::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
