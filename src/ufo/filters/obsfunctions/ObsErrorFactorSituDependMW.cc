/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSituDependMW.h"

#include <cmath>

#include <algorithm>
#include <memory>
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

static ObsFunctionMaker<ObsErrorFactorSituDependMW>
       makerObsErrorFactorSituDependMW_("ObsErrorFactorSituDependMW");

// -----------------------------------------------------------------------------

ObsErrorFactorSituDependMW::ObsErrorFactorSituDependMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);
  ASSERT(inst == "amsua" || inst == "atms");

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Include list of required data from ObsSpace
  invars_ += Variable("ObsValue/brightnessTemperature", channels_);
  invars_ += Variable("ObsError/brightnessTemperature", channels_);
  invars_ += Variable(hofxgrp+"/brightnessTemperature", channels_);
  invars_ += Variable(errgrp+"/brightnessTemperature", channels_);
  invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
  invars_ += Variable("GeoVaLs/surface_wind_speed");

  // Include required variables from ObsFunction
  const Variable &clwobs = options_.clwobsFunction.value();
  invars_ += clwobs;

  const Variable &clwbkg = options_.clwbkgFunction.value();
  invars_ += clwbkg;

  const Variable &scatobs = options_.scatobsFunction.value();
  invars_ += scatobs;

  const Variable &clwmatchidx = options_.clwmatchidxFunction.value();
  invars_ += clwmatchidx;

  if (options_.obserrFunction.value() != boost::none) {
    const boost::optional<Variable> &obserrvar = options_.obserrFunction.value();
    invars_ += *obserrvar;
  }
}

// -----------------------------------------------------------------------------

ObsErrorFactorSituDependMW::~ObsErrorFactorSituDependMW() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorSituDependMW::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get test groups from options
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Get clear-sky observation error from options
  const std::vector<float> &obserr_clr = options_.obserrClearSky.value();

  // Get CLW retrieval based on observation from ObsFunction
  const Variable &clwobsvar = options_.clwobsFunction.value();
  ioda::ObsDataVector<float> clwobs(in.obsspace(), clwobsvar.toOopsObsVariables());
  in.get(clwobsvar, clwobs);

  // Get CLW retrieval based on simulated observation from ObsFunction
  const Variable &clwbkgvar = options_.clwbkgFunction.value();
  ioda::ObsDataVector<float> clwbkg(in.obsspace(), clwbkgvar.toOopsObsVariables());
  in.get(clwbkgvar, clwbkg);

  // Get Scattering Index retrieval based on observation from ObsFunction
  const Variable &scatobsvar = options_.scatobsFunction.value();
  ioda::ObsDataVector<float> scatobs(in.obsspace(), scatobsvar.toOopsObsVariables());
  in.get(scatobsvar, scatobs);

  // Get CLW match index from ObsFunction
  const Variable &clwmatchidxvar = options_.clwmatchidxFunction.value();
  ioda::ObsDataVector<float> clwmatchidx(in.obsspace(), clwmatchidxvar.toOopsObsVariables());
  in.get(clwmatchidxvar, clwmatchidx);

  // Get ObsErrorData (obs error from previous QC step) and convert to inverse of error variance
  std::vector<std::vector<float>> varinv(nchans, std::vector<float>(nlocs));
  std::vector<float> obserrdata;
  std::vector<int> qcflagdata;
  const float missing = util::missingValue<float>();
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);
    in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
    for (size_t iloc = 0; iloc < nlocs; iloc++) {
      if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                           : qcflagdata[iloc] = 0;
      (qcflagdata[iloc] == 0) ? (varinv[ichan][iloc] = 1.0 / pow(obserrdata[iloc], 2))
                              : (varinv[ichan][iloc] = 0.0);
    }
  }

  // Calculate bias corrected innovation: ObsValue - HofX (HofX includes bias correction)
  std::vector<std::vector<float>> btobs(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> hofx(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> innov(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("ObsValue/brightnessTemperature", channels_)[ichan], btobs[ichan]);
    in.get(Variable(hofxgrp+"/brightnessTemperature", channels_)[ichan], hofx[ichan]);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innov[ichan][iloc] = btobs[ichan][iloc] - hofx[ichan][iloc];
    }
  }

  // Get variables from GeoVaLs
  // Get surface wind speed
  std::vector<float> surface_wind_speed(nlocs);
  in.get(Variable("GeoVaLs/surface_wind_speed"), surface_wind_speed);

  // Load area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  // Set channel number
  int ich238, ich314, ich503, ich528, ich536, ich544, ich549, ich890;
  if (inst == "amsua") {
    ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 4, ich536 = 5;
    ich544 = 6, ich549 = 7, ich890 = 15;
  } else if (inst == "atms") {
    ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 5, ich536 = 6;
    ich544 = 7, ich549 = 8, ich890 = 16;
  }

  // Get Original Observation Error from ObsFunction
  std::unique_ptr<ioda::ObsDataVector<float>> obserr0;
  if (options_.obserrFunction.value() != boost::none) {
    const boost::optional<Variable> &obserrvar = options_.obserrFunction.value();
    obserr0.reset(new ioda::ObsDataVector<float>(in.obsspace(),
                 (*obserrvar).toOopsObsVariables()));
    in.get(*obserrvar, *obserr0);
    // Calculate error factors (error_factors) for each channel
    // Loop through locations
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) out[ichan][iloc] = 1.0;
      if (water_frac[iloc] >= 0.99) {
        float icol = 1.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) icol = icol * clwmatchidx[ichan][iloc];
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          size_t channel = ichan + 1;
          if (varinv[ichan][iloc] > 0.0 && (channel <= ich536 || channel >= ich890)) {
            float term = (1.0 - icol) * std::abs(innov[ichan][iloc]);
            term = term + std::min(0.002 * pow(surface_wind_speed[iloc], 2) *
                          (*obserr0)[ichan][iloc], 0.5 * (*obserr0)[ichan][iloc]);
            float clwtmp = std::min(std::abs((clwobs[0][iloc] - clwbkg[0][iloc])), 1.f);
            term = term + std::min(13.0 * clwtmp * (*obserr0)[ichan][iloc], 3.5 *
                          (*obserr0)[ichan][iloc]);
            if (scatobs[0][iloc] > 9.0) {
              term = term + std::min(1.5 * (scatobs[0][iloc] - 9.0) * (*obserr0)[ichan][iloc],
                                     2.5 * (*obserr0)[ichan][iloc]);
            }
            term = pow(term, 2.0);
            out[ichan][iloc] = 1.0 / (1.0 + varinv[ichan][iloc] * term);
            out[ichan][iloc] = sqrt(1.0 / out[ichan][iloc]);
          }
        }
      }
    }
  } else {
    std::vector<std::vector<float>> obserr0(nchans, std::vector<float>(nlocs));
    // Get Original Observation Error (if ObsError is filled up with missing values, replace it)
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      in.get(Variable("ObsError/brightnessTemperature", channels_)[ichan], obserr0[ichan]);
      for (size_t iloc = 0; iloc < nlocs; iloc++) {
        if (obserr0[ichan][iloc] == missing) obserr0[ichan][iloc] = obserr_clr[ichan];
      }
    }
    // Calculate error factors (error_factors) for each channel
    // Loop through locations
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) out[ichan][iloc] = 1.0;
      if (water_frac[iloc] >= 0.99) {
        float icol = 1.0;
        for (size_t ichan = 0; ichan < nchans; ++ichan) icol = icol * clwmatchidx[ichan][iloc];
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          size_t channel = ichan + 1;
          if (varinv[ichan][iloc] > 0.0 && (channel <= ich536 || channel >= ich890)) {
            float term = (1.0 - icol) * std::abs(innov[ichan][iloc]);
            term = term + std::min(0.002 * pow(surface_wind_speed[iloc], 2) *
                          obserr0[ichan][iloc], 0.5 * obserr0[ichan][iloc]);
            float clwtmp = std::min(std::abs((clwobs[0][iloc] - clwbkg[0][iloc])), 1.f);
            term = term + std::min(13.0 * clwtmp * obserr0[ichan][iloc], 3.5 *
                          obserr0[ichan][iloc]);
            if (scatobs[0][iloc] > 9.0) {
              term = term + std::min(1.5 * (scatobs[0][iloc] - 9.0) * obserr0[ichan][iloc],
                                     2.5 * obserr0[ichan][iloc]);
            }
            term = pow(term, 2.0);
            out[ichan][iloc] = 1.0 / (1.0 + varinv[ichan][iloc] * term);
            out[ichan][iloc] = sqrt(1.0 / out[ichan][iloc]);
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSituDependMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
