/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/HydrometeorCheckAMSUA.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/filters/obsfunctions/ObsErrorModelRamp.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"


namespace ufo {

static ObsFunctionMaker<HydrometeorCheckAMSUA> makerHydrometeorCheckAMSUA_("HydrometeorCheckAMSUA");

// -----------------------------------------------------------------------------

HydrometeorCheckAMSUA::HydrometeorCheckAMSUA(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &biastermgrp_ = options_.testBiasTerm.value();
  const std::string &biasgrp_ = options_.testBias.value();
  const std::string &hofxgrp_ = options_.testHofX.value();

  // Include required variables from ObsDiag
  invars_ += Variable("brightness_temperature_jacobian_surface_emissivity@ObsDiag", channels_);
  invars_ += Variable("brightness_temperature_assuming_clear_sky@ObsDiag", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@"+biasgrp_, channels_);
  invars_ += Variable("brightness_temperature@"+hofxgrp_, channels_);
  invars_ += Variable("constant@"+biastermgrp_, channels_);
  invars_ += Variable("scan_angle_4th_order@"+biastermgrp_, channels_);
  invars_ += Variable("scan_angle_3rd_order@"+biastermgrp_, channels_);
  invars_ += Variable("scan_angle_2nd_order@"+biastermgrp_, channels_);
  invars_ += Variable("scan_angle_1st_order@"+biastermgrp_, channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");

  // Include list of required data from ObsFunction
  const Variable &obserrfunc = options_.obserrFunction.value();
  invars_ += obserrfunc;

  const Variable &clwretfunc = options_.clwretFunction.value();
  invars_ += clwretfunc;
}

// -----------------------------------------------------------------------------

HydrometeorCheckAMSUA::~HydrometeorCheckAMSUA() {}

// -----------------------------------------------------------------------------

void HydrometeorCheckAMSUA::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get test groups from options
  const std::string &biastermgrp_ = options_.testBiasTerm.value();
  const std::string &biasgrp_ = options_.testBias.value();
  const std::string &hofxgrp_ = options_.testHofX.value();

  // Get clear-sky observation error from options
  const std::vector<float> &obserr0 = options_.obserrClearSky.value();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Get surface temperature jacobian
  std::vector<std::vector<float>> dbtde(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature_jacobian_surface_emissivity@ObsDiag", channels_)[ichan],
           dbtde[ichan]);
  }

  // Get Observation
  std::vector<std::vector<float>> btobs(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ichan], btobs[ichan]);
  }

  // Get Observation Bias
  std::vector<std::vector<float>> bias(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+biasgrp_, channels_)[ichan], bias[ichan]);
  }

  // Get HofX
  std::vector<std::vector<float>> hofx(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+hofxgrp_, channels_)[ichan], hofx[ichan]);
  }

  // Get HofX for clear-sky simulation
  std::vector<std::vector<float>> hofxclr(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature_assuming_clear_sky@ObsDiag", channels_)[ichan],
           hofxclr[ichan]);
  }

  // Get ObsBiasTerm: constant term
  std::vector<std::vector<float>> bias_const(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("constant@"+biastermgrp_, channels_)[ichan], bias_const[ichan]);
  }

  // Get ObsBiasTerm: scan angle terms
  size_t nangs = 4;
  std::vector<float> values(nlocs);
  std::vector<std::string> scanterms(nangs);
  std::vector<std::vector<float>> bias_scanang(nchans, std::vector<float>(nlocs));
  scanterms[0] = "scan_angle_4th_order@"+biastermgrp_;
  scanterms[1] = "scan_angle_3rd_order@"+biastermgrp_;
  scanterms[2] = "scan_angle_2nd_order@"+biastermgrp_;
  scanterms[3] = "scan_angle_1st_order@"+biastermgrp_;
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iang = 0; iang < nangs; ++iang) {
      in.get(Variable(scanterms[iang], channels_)[ichan], values);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        bias_scanang[ichan][iloc] = bias_scanang[ichan][iloc] + values[iloc];
      }
    }
  }

  // Calculate bias-corrected innovation: Observation - HofX - bias
  std::vector<std::vector<float>> innov(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innov[ichan][iloc] = btobs[ichan][iloc] - hofx[ichan][iloc] - bias[ichan][iloc];
    }
  }

  // Get all-sky observation error from ObsFunction
  const Variable &obserrvar = options_.obserrFunction.value();
  ioda::ObsDataVector<float> obserr(in.obsspace(), obserrvar.toOopsVariables());
  in.get(obserrvar, obserr);

  // Get CLW retrieval based on observation from ObsFunction
  const Variable &clwretvar = options_.clwretFunction.value();
  ioda::ObsDataVector<float> clwobs(in.obsspace(), clwretvar.toOopsVariables());
  in.get(clwretvar, clwobs);

  // Set channel index
  int ich238 = 0, ich314 = 1, ich503 = 2, ich528 = 3, ich536 = 4;
  int ich544 = 5, ich549 = 6, ich890 = 14;

  // Set parameters
  float w1f6 = 1.0/10.0, w2f6 = 1.0/0.80;
  float w1f4 = 1.0/0.30, w2f4 = 1.0/1.80;

  std::vector<std::vector<int>> affected_channels(nchans, std::vector<int>(nlocs));
  // Loop over locations
  // Combined cloud-precipitation-surface checks
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // Initialization
    for (size_t ich = 0; ich < nchans; ++ich) {
      affected_channels[ich][iloc] = 0;
    }

    // Calculate cloud effect
    float cldeff_obs536 = 0.0;
    if (water_frac[iloc] > 0.99) {
      cldeff_obs536 = btobs[ich536][iloc] - hofxclr[ich536][iloc] - bias[ich536][iloc];
    }

    // Calculate scattering effect
    std::vector<float> factch4(nlocs);
    std::vector<float> factch6(nlocs);
    float btobsbc238 = btobs[ich238][iloc] - bias_const[ich238][iloc]
                                           - bias_scanang[ich238][iloc];
    float clwx = 0.6;
    float dsval = 0.8;
    if (water_frac[iloc] > 0.99) {
      clwx = 0.0;
      dsval = ((2.410 - 0.0098 * btobsbc238) * innov[ich238][iloc] +
                0.454 * innov[ich314][iloc] - innov[ich890][iloc]) * w1f6;
      dsval = std::max(static_cast<float>(0.0), dsval);
    }
    factch4[iloc] = pow(clwx, 2) + pow(innov[ich528][iloc] * w2f4, 2);
    factch6[iloc] = pow(dsval, 2) + pow(innov[ich544][iloc] * w2f6, 2);

    // Window channel sanity check
    // If any of the window channels is bad, skip all window channels
    // List of surface sensitivity channels
    std::vector<float> OmFs{std::abs(innov[ich238][iloc]), std::abs(innov[ich314][iloc]),
                            std::abs(innov[ich503][iloc]), std::abs(innov[ich528][iloc]),
                            std::abs(innov[ich536][iloc]), std::abs(innov[ich544][iloc]),
                            std::abs(innov[ich890][iloc])};
    bool result = false;
    result = any_of(OmFs.begin(), OmFs.end(), [](float x){return x > 200.0;});

    if (result) {
      for (size_t ich = ich238; ich <= ich544; ++ich) {
        affected_channels[ich][iloc] = 1;
      }
      affected_channels[ich890][iloc] = 1;
    } else {
      // Hydrometeor check over water surface
      if (water_frac[iloc] > 0.99) {
        // Cloud water retrieval sanity check
        if (clwobs[0][iloc] > 999.0) {
          for (size_t ich = ich238; ich <= ich544; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          affected_channels[ich890][iloc] = 1;
        }

        // Precipitation check (factch6)
        if (factch6[iloc] >= 1.0) {
          for (size_t ich = ich238; ich <= ich544; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          affected_channels[ich890][iloc] = 1;
        // Scattering check (ch5 cloud effect)
        } else if (cldeff_obs536 < -0.5) {
          for (size_t ich = ich238; ich <= ich544; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          affected_channels[ich890][iloc] = 1;
        // Sensitivity of BT to the surface emissivity check
        } else {
          float thrd238 = 0.025, thrd314 = 0.015, thrd503 = 0.030, thrd890 = 0.030;
          float de238 = 0.0, de314 = 0.0, de503 = 0.0, de890 = 0.0;
          float dbtde238 = dbtde[ich238][iloc];
          float dbtde314 = dbtde[ich314][iloc];
          float dbtde503 = dbtde[ich503][iloc];
          float dbtde890 = dbtde[ich890][iloc];
          if (dbtde238 != 0.0) de238 = std::abs(innov[ich238][iloc]) / dbtde238 *
                                       (obserr0[ich238] / obserr[ich238][iloc]) *
                                       (1.0 - std::max(1.0, 10.0*clwobs[0][iloc]));
          if (dbtde314 != 0.0) de314 = std::abs(innov[ich314][iloc]) / dbtde314 *
                                       (obserr0[ich314] / obserr[ich314][iloc]) *
                                       (1.0 - std::max(1.0, 10.0*clwobs[0][iloc]));
          if (dbtde503 != 0.0) de503 = std::abs(innov[ich503][iloc]) / dbtde503 *
                                       (obserr0[ich503] / obserr[ich503][iloc]) *
                                       (1.0 - std::max(1.0, 10.0*clwobs[0][iloc]));
          if (dbtde890 != 0.0) de890 = std::abs(innov[ich890][iloc]) / dbtde890 *
                                       (obserr0[ich890] / obserr[ich890][iloc]) *
                                       (1.0 - std::max(1.0, 10.0*clwobs[0][iloc]));
          bool qcemiss = false;
          qcemiss = de238 > thrd238 || de314 > thrd314 || de503 > thrd503 || de890 > thrd890;
          if (qcemiss) {
            for (size_t ich = ich238; ich <= ich536; ++ich) {
              affected_channels[ich][iloc] = 1;
            }
            affected_channels[ich890][iloc] = 1;
          }
        }
      } else {
        // Hydrometeor check over non-water (land/sea ice/snow) surface
        // Precipitation check (factch6)
        if (factch6[iloc] >= 1.0) {
          for (size_t ich = ich238; ich <= ich544; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          affected_channels[ich890][iloc] = 1;
        // Thick cloud check (factch4)
        } else if (factch4[iloc] > 0.5) {
          for (size_t ich = ich238; ich <= ich536; ++ich) {
            affected_channels[ich][iloc] = 1;
          }
          affected_channels[ich890][iloc] = 1;
        // Sensitivity of BT to the surface emissivity check
        } else {
          float thrd238 = 0.020, thrd314 = 0.015, thrd503 = 0.035, thrd890 = 0.015;
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
            affected_channels[ich890][iloc] = 1;
          }
        }
      // surface type
      }
    // window channel sanity check
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

const ufo::Variables & HydrometeorCheckAMSUA::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
