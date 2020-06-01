/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSituDependMW.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/CLWMatchIndexMW.h"
#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/filters/obsfunctions/SCATRetMW.h"
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

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();
  const std::string &biasgrp = options_.testBias.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@"+hofxgrp, channels_);
  invars_ += Variable("brightness_temperature@"+biasgrp, channels_);
  invars_ += Variable("brightness_temperature@"+errgrp, channels_);
  invars_ += Variable("brightness_temperature@"+flaggrp, channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("surface_wind_speed@GeoVaLs");

  // Include required variables from ObsFunction
  const Variable &clwobs = options_.clwobsFunction.value();
  invars_ += clwobs;

  const Variable &clwbkg = options_.clwbkgFunction.value();
  invars_ += clwbkg;

  const Variable &scatobs = options_.scatobsFunction.value();
  invars_ += scatobs;

  const Variable &clwmatchidx = options_.clwmatchidxFunction.value();
  invars_ += clwmatchidx;
}

// -----------------------------------------------------------------------------

ObsErrorFactorSituDependMW::~ObsErrorFactorSituDependMW() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorSituDependMW::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get test groups from options
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &biasgrp = options_.testBias.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Get clear-sky observation error from options
  const std::vector<float> &obserr_clr = options_.obserrClearSky.value();

  // Get CLW retrieval based on observation from ObsFunction
  const Variable &clwobsvar = options_.clwobsFunction.value();
  ioda::ObsDataVector<float> clwobs(in.obsspace(), clwobsvar.toOopsVariables());
  in.get(clwobsvar, clwobs);

  // Get CLW retrieval based on simulated observation from ObsFunction
  const Variable &clwbkgvar = options_.clwbkgFunction.value();
  ioda::ObsDataVector<float> clwbkg(in.obsspace(), clwbkgvar.toOopsVariables());
  in.get(clwbkgvar, clwbkg);

  // Get Scattering Index retrieval based on observation from ObsFunction
  const Variable &scatobsvar = options_.scatobsFunction.value();
  ioda::ObsDataVector<float> scatobs(in.obsspace(), scatobsvar.toOopsVariables());
  in.get(scatobsvar, scatobs);

  // Get CLW match index from ObsFunction
  const Variable &clwmatchidxvar = options_.clwmatchidxFunction.value();
  ioda::ObsDataVector<float> clwmatchidx(in.obsspace(), clwmatchidxvar.toOopsVariables());
  in.get(clwmatchidxvar, clwmatchidx);

  // Get Original Observation Error
  std::vector<std::vector<float>> obserr0(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsError", channels_)[ichan], obserr0[ichan]);
  }

  // Get ObsErrorData (obs error from previous QC step) and convert to inverse of error variance
  std::vector<std::vector<float>> varinv(nchans, std::vector<float>(nlocs));
  std::vector<float> obserrdata;
  std::vector<int> qcflagdata;
  const float missing = util::missingValue(missing);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+flaggrp, channels_)[ichan], qcflagdata);
    in.get(Variable("brightness_temperature@"+errgrp, channels_)[ichan], obserrdata);
    for (size_t iloc = 0; iloc < nlocs; iloc++) {
      if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                           : qcflagdata[iloc] = 0;
      (qcflagdata[iloc] == 0) ? (varinv[ichan][iloc] = 1.0 / pow(obserrdata[iloc], 2))
                              : (varinv[ichan][iloc] = 0.0);
    }
  }

  // Calculate bias corrected innovation: ObsValue - HofX - ObsBias
  std::vector<std::vector<float>> btobs(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> hofx(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> bias(nchans, std::vector<float>(nlocs));
  std::vector<std::vector<float>> innov(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ichan], btobs[ichan]);
    in.get(Variable("brightness_temperature@"+hofxgrp, channels_)[ichan], hofx[ichan]);
    in.get(Variable("brightness_temperature@"+biasgrp, channels_)[ichan], bias[ichan]);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innov[ichan][iloc] = btobs[ichan][iloc] - hofx[ichan][iloc] - bias[ichan][iloc];
    }
  }

  // Get variables from GeoVaLs
  // Get surface wind speed
  std::vector<float> surface_wind_speed(nlocs);
  in.get(Variable("surface_wind_speed@GeoVaLs"), surface_wind_speed);

  // Load area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Calculate error factors (error_factors) for each channel
  // Set channel number
  int ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 4, ich536 = 5;
  int ich544 = 6, ich549 = 7, ich890 = 15;
  // Loop through locations
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    for (size_t ichan = 0; ichan < nchans; ++ichan) out[ichan][iloc] = 1.0;
    if (water_frac[iloc] >= 0.99) {
      float icol = 1.0;
      for (size_t ichan = 0; ichan < nchans; ++ichan) icol = icol * clwmatchidx[ichan][iloc];
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        size_t channel = ichan + 1;
        if (varinv[ichan][iloc] > 0.0 && (channel <= ich536 || channel == ich890)) {
          float term = (1.0 - icol) * std::abs(innov[ichan][iloc]);
          term = term + std::min(0.002 * pow(surface_wind_speed[iloc], 2) * obserr0[ichan][iloc],
                                 0.5 * obserr0[ichan][iloc]);
          float clwtmp = std::min(std::abs((clwobs[0][iloc] - clwbkg[0][iloc])), 1.f);
          term = term + std::min(13.0 * clwtmp * obserr0[ichan][iloc], 3.5 * obserr0[ichan][iloc]);
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

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSituDependMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
