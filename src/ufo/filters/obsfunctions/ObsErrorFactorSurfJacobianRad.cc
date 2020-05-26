/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSurfJacobianRad.h"

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

static ObsFunctionMaker<ObsErrorFactorSurfJacobianRad>
       makerObsErrorFactorSurfJacobianRad_("ObsErrorFactorSurfJacobianRad");

// -----------------------------------------------------------------------------

ObsErrorFactorSurfJacobianRad::ObsErrorFactorSurfJacobianRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();

  // Include required variables from ObsDiag
  invars_ += Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag", channels_);
  invars_ += Variable("brightness_temperature_jacobian_surface_emissivity@ObsDiag", channels_);

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@"+errgrp, channels_);
  invars_ += Variable("brightness_temperature@"+flaggrp, channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("land_area_fraction@GeoVaLs");
  invars_ += Variable("ice_area_fraction@GeoVaLs");
  invars_ += Variable("surface_snow_area_fraction@GeoVaLs");
}

// -----------------------------------------------------------------------------

ObsErrorFactorSurfJacobianRad::~ObsErrorFactorSurfJacobianRad() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorSurfJacobianRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get options from configuration
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Load area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  std::vector<float> land_frac(nlocs);
  std::vector<float> ice_frac(nlocs);
  std::vector<float> snow_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);
  in.get(Variable("land_area_fraction@GeoVaLs"), land_frac);
  in.get(Variable("ice_area_fraction@GeoVaLs"), ice_frac);
  in.get(Variable("surface_snow_area_fraction@GeoVaLs"), snow_frac);

  // Get observation tuning parameters over sea/land/oce/snow/mixed from options
  const std::vector<float> &demisf_in = options_.obserrScaleFactorEsfc.value();
  const std::vector<float> &dtempf_in = options_.obserrScaleFactorTsfc.value();

  // Setup weight for each surface type
  std::vector<float> demisf(nlocs);
  std::vector<float> dtempf(nlocs);

  // Determine surface type and weight for current obs
  for (size_t iloc = 0; iloc < nlocs; iloc++) {
    bool sea = water_frac[iloc] >= 0.99;
    bool land = land_frac[iloc] >= 0.99;
    bool ice = ice_frac[iloc] >= 0.99;
    bool snow = snow_frac[iloc] >= 0.99;
    bool mixed = (!sea && !land && !ice && !snow);
    if (sea) {
      demisf[iloc] = demisf_in[0];
      dtempf[iloc] = dtempf_in[0];
    }
    if (land) {
      demisf[iloc] = demisf_in[1];
      dtempf[iloc] = dtempf_in[1];
    }
    if (ice) {
      demisf[iloc] = demisf_in[2];
      dtempf[iloc] = dtempf_in[2];
    }
    if (snow) {
      demisf[iloc] = demisf_in[3];
      dtempf[iloc] = dtempf_in[3];
    }
    if (mixed) {
      demisf[iloc] = demisf_in[4];
      dtempf[iloc] = dtempf_in[4];
    }
  }

  // Calculate error factors for each channel
  std::vector<int> qcflagdata;
  std::vector<float> obserrdata;
  std::vector<float> dbtdts(nlocs);
  std::vector<float> dbtdes(nlocs);
  float varinv = 0.0;
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const float missing = util::missingValue(missing);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag",
                     channels_)[ichan], dbtdts);
    in.get(Variable("brightness_temperature_jacobian_surface_emissivity@ObsDiag",
                     channels_)[ichan], dbtdes);
    in.get(Variable("brightness_temperature@"+errgrp, channels_)[ichan], obserrdata);
    in.get(Variable("brightness_temperature@"+flaggrp, channels_)[ichan], qcflagdata);

    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                           : qcflagdata[iloc] = 0;
      (qcflagdata[iloc] == 0) ? (varinv = 1.0 / pow(obserrdata[iloc], 2)) : (varinv = 0.0);
      out[ichan][iloc] = 1.0;

      if (varinv > 0.0) {
        float vaux = demisf[iloc] * std::fabs(dbtdes[iloc]) +
               dtempf[iloc] * std::fabs(dbtdts[iloc]);
        float term = pow(vaux, 2);
        if (term > 0.0) {
          out[ichan][iloc] = sqrt(1.0 / (1.0 / (1.0 + varinv * term)));
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSurfJacobianRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
