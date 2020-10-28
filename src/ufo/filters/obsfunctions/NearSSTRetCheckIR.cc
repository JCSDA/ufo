/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/NearSSTRetCheckIR.h"

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
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<NearSSTRetCheckIR> makerNearSSTRetCheckIR_("NearSSTRetCheckIR");

// -----------------------------------------------------------------------------

NearSSTRetCheckIR::NearSSTRetCheckIR(const eckit::LocalConfiguration & conf)
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
  const std::string &biasgrp = options_.testBias.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Include required variables from ObsDiag
  invars_ += Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag", channels_);
  invars_ += Variable("brightness_temperature_jacobian_air_temperature@ObsDiag", channels_);
  invars_ += Variable("brightness_temperature_jacobian_humidity_mixing_ratio@ObsDiag", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@"+flaggrp, channels_);
  invars_ += Variable("brightness_temperature@"+errgrp, channels_);
  invars_ += Variable("brightness_temperature@"+hofxgrp, channels_);
  invars_ += Variable("brightness_temperature@"+biasgrp, channels_);
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@ObsError", channels_);
  invars_ += Variable("sensor_band_central_radiation_wavenumber@VarMetaData");

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("surface_temperature_where_sea@GeoVaLs");
}

// -----------------------------------------------------------------------------

NearSSTRetCheckIR::~NearSSTRetCheckIR() {}

// -----------------------------------------------------------------------------

void NearSSTRetCheckIR::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get channel usage information from options
  std::vector<int> use_flag = options_.useflagChannel.value();

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("brightness_temperature_jacobian_air_temperature@ObsDiag",
                                    channels_)[0]);

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();
  const std::string &biasgrp = options_.testBias.value();
  const std::string &hofxgrp = options_.testHofX.value();

  // Setup vectors to get 2D variables
  std::vector<float> values(nlocs);

  // Get variables from ObsDiag
  // Get surface temperature jacobian
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

  // Get moisture jacobian
  std::vector<std::vector<std::vector<float>>>
       dbtdq(nchans, std::vector<std::vector<float>>(nlevs, std::vector<float>(nlocs)));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      int level = nlevs - ilev;
      in.get(Variable("brightness_temperature_jacobian_humidity_mixing_ratio@ObsDiag",
                       channels_)[ichan], level, dbtdq[ichan][ilev]);
    }
  }

  // Get variables from ObsSpace
  // Get sensor band central radiation wavenumber
  std::vector<float> wavenumber(nchans);
  in.get(Variable("sensor_band_central_radiation_wavenumber@VarMetaData"), wavenumber);

  // Get effective observation error and convert it to inverse of the error variance
  const float missing = util::missingValue(missing);
  std::vector<int> qcflag(nlocs, 0);
  std::vector<std::vector<float>> varinv(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+errgrp, channels_)[ichan], values);
    in.get(Variable("brightness_temperature@"+flaggrp, channels_)[ichan], qcflag);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (flaggrp == "PreQC") values[iloc] == missing ? qcflag[iloc] = 100 : qcflag[iloc] = 0;
      (qcflag[iloc] == 0) ? (varinv[ichan][iloc] = 1.0 / pow(values[iloc], 2))
                          : (varinv[ichan][iloc] = 0.0);
    }
  }

  // Get bias corrected innovation (tbobs - hofx - bias)
  std::vector<std::vector<float>> innovation(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ichan], innovation[ichan]);
    in.get(Variable("brightness_temperature@"+hofxgrp, channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innovation[ichan][iloc] = innovation[ichan][iloc] - values[iloc];
    }
    in.get(Variable("brightness_temperature@"+biasgrp, channels_)[ichan], values);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      innovation[ichan][iloc] = innovation[ichan][iloc] - values[iloc];
    }
  }

  // Get original observation error (uninflated)
  std::vector<std::vector<float>> obserr(nchans, std::vector<float>(nlocs));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsError", channels_)[ichan], obserr[ichan]);
  }

  // Get variables from GeoVaLS
  // Get solar zenith angle
  std::vector<float> solza(nlocs);
  in.get(Variable("solar_zenith_angle@MetaData"), solza);

  // Get water temperature
  std::vector<float> tzbgr(nlocs);
  in.get(Variable("surface_temperature_where_sea@GeoVaLs"), tzbgr);

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Retrieved NSST increment (dtz) and average of surface temeprature jacobian (ts_ave)
  // Setup constants
  // tschk: threshold for surface temperature jacobian
  // tzchk: threshold for SST temperature at obs location
  const float tschk = 0.2, tzchk = 0.85;
  const float e_ts = 0.5, e_ta = 1.f, e_qa = 0.85;
  const float ttp = Constants::ttp;

  // Loop through locations
  // TODO(EL): review the use of irday with EMC
  std::vector<int> irday(nchans, 1);
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    bool sea = water_frac[iloc] >= 0.99;
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      out[ichan][iloc] = 0.0;
      if (water_frac[iloc] > 0.0 && solza[iloc] <= 89.0 && wavenumber[ichan] > 2400.0) {
        irday[ichan] = 0;
      }
    }
    int icount = 0;
    float ws, wa, wq;
    float a11 = 0.0, a22 = 0.0, a33 = 0.0, a12 = 0.0, a13 = 0.0, a23 = 0.0;
    float c1x = 0.0, c2x = 0.0, c3x = 0.0;
    float dtz = -999.0;
    float ts_ave = 0.0;
    std::vector<float> tb_ta(nchans);
    std::vector<float> tb_qa(nchans);
    if (sea) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        if (use_flag[ichan] >= 1 && varinv[ichan][iloc] > 0.0 && irday[ichan] == 1
                                 && dbtdts[ichan][iloc] >= tschk) {
          tb_ta[ichan] = dbtdt[ichan][nlevs-1][iloc];
          tb_qa[ichan] = dbtdq[ichan][nlevs-1][iloc];
          for (size_t ilev = 0; ilev < nlevs-1; ++ilev) {
            tb_ta[ichan] = tb_ta[ichan] + dbtdt[ichan][ilev][iloc];
            tb_qa[ichan] = tb_qa[ichan] + dbtdq[ichan][ilev][iloc];
          }
        }
      }
      ws = 1.0 / pow(e_ts, 2);
      wa = 1.0 / pow(e_ta, 2);
      wq = 1.0 / pow(e_qa * (std::max(((tzbgr[iloc]) - ttp)
                          * 0.03, 0.0) + 0.1), 2);
      a11 = ws;
      a22 = wa;
      a33 = wq;

      // Get coefficients for linear equations
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        if (use_flag[ichan] >= 1 && varinv[ichan][iloc] > 0.0 && irday[ichan] == 1
                                 && dbtdts[ichan][iloc] >= tschk) {
          icount = icount + 1;
          ts_ave = ts_ave + dbtdts[ichan][iloc];
          float w_rad = pow((1.0 / obserr[ichan][iloc]), 2);
          a11 = a11 + w_rad * dbtdts[ichan][iloc] * dbtdts[ichan][iloc];
          a12 = a12 + w_rad * dbtdts[ichan][iloc] * tb_ta[ichan];
          a13 = a13 + w_rad * dbtdts[ichan][iloc] * tb_qa[ichan];
          a22 = a22 + w_rad * tb_ta[ichan] * tb_ta[ichan];
          float varrad = w_rad * innovation[ichan][iloc];
          c1x = c1x + varrad * dbtdts[ichan][iloc];
          c2x = c2x + varrad * tb_ta[ichan];
          c3x = c3x + varrad * tb_qa[ichan];
        }
      }

      // Solve linear equations with three unknowns (dtz, dta, dqa)
      // only dtz is solved since other two are not useful here
      float delt1 = 0.0, delt = 1.0;
      if ( icount >= 1 ) {
        delt = a11 * (a22 * a33 - a23 * a23) +
               a12 * (a13 * a23 - a12 * a33) +
               a13 * (a12 * a23 - a13 * a22);
        delt1 = c1x * (a22 * a33 - a23 * a23) +
                c2x * (a13 * a23 - a12 * a33) +
                c3x * (a12 * a23 - a13 * a22);
        dtz = delt1 / delt;
        ts_ave = ts_ave / static_cast<float>(icount);
      }
    // sea
    }

    if (dtz != -999.0) {
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        if (use_flag[ichan] >= 1 && varinv[ichan][iloc] > 0.0 && dbtdts[ichan][iloc] >= tschk) {
          float xindx = pow((dbtdts[ichan][iloc] - ts_ave) / (1.0 - ts_ave), 3);
          float tzchks = tzchk * pow(0.5, xindx);
          if (std::fabs(dtz) > tzchks) out[ichan][iloc] = 1;
        }
      }
    }
  // end of location loop
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & NearSSTRetCheckIR::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
