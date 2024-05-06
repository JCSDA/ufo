/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SurfTypeCheckRad.h"

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

static ObsFunctionMaker<SurfTypeCheckRad>
       makerSurfTypeCheckRad_("SurfTypeCheckRad");

// -----------------------------------------------------------------------------

SurfTypeCheckRad::SurfTypeCheckRad(const eckit::LocalConfiguration & conf)
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

  // Include list of required data from ObsSpace
  invars_ += Variable(errgrp+"/brightnessTemperature", channels_);
  invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
  invars_ += Variable("GeoVaLs/land_area_fraction");
  invars_ += Variable("GeoVaLs/ice_area_fraction");
  invars_ += Variable("GeoVaLs/surface_snow_area_fraction");
}

// -----------------------------------------------------------------------------

SurfTypeCheckRad::~SurfTypeCheckRad() {}

// -----------------------------------------------------------------------------

void SurfTypeCheckRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get channel use flags from options
  std::vector<int> use_flag = options_.useflagChannel.value();
  std::vector<int> use_flag_clddet = options_.useflagCloudDetect.value();

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

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
  size_t iwater_det = 0;
  size_t iland_det = 0;
  size_t isnow_det = 0;
  size_t iice_det = 0;
  size_t imix_det = 0;
  size_t ndim = 6;
  int dec = 0;
  int bindec = 0;
  std::vector<int> qcflagdata(nlocs);
  std::vector<float> obserrdata(nlocs);
  std::vector<int> bin(ndim, 0);
  float varinv = 0.0;
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const float missing = util::missingValue<float>();

  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
    in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);

    for (size_t iloc = 0; iloc < nlocs; ++iloc) out[ichan][iloc] = 0.0;

    if (use_flag[ichan] >= 1 && use_flag_clddet[ichan] >= 2) {
       ASSERT((use_flag_clddet[ichan] - pow(2, ndim)) < 0);

       dec = use_flag_clddet[ichan];
       for (size_t i = ndim; i >= 1; --i) {
          bindec = pow(2, i-1);
          if (dec >= bindec) {
            bin[i-1] = 1;
            dec = dec - bindec;
          } else {
            bin[i-1] = 0;
          }
       }
       iland_det = bin[1];
       isnow_det = bin[2];
       imix_det = bin[3];
       iice_det = bin[4];
       iwater_det = bin[5];

       for (size_t iloc = 0; iloc < nlocs; ++iloc) {
         sea = water_frac[iloc] >= 0.99;
         land = land_frac[iloc] >= 0.99;
         ice = ice_frac[iloc] >= 0.99;
         snow = snow_frac[iloc] >= 0.99;
         mixed = (!sea && !land && !ice && !snow);
         if (flaggrp == "PreQC")
                  obserrdata[iloc] == missing ? qcflagdata[iloc] = 100 : qcflagdata[iloc] = 0;
         (qcflagdata[iloc] == 0) ? (varinv = 1.0 / pow(obserrdata[iloc], 2)) : (varinv = 0.0);
         if (varinv > 0.0) {
           if (sea && iwater_det > 0)  out[ichan][iloc] = 1.0;
           else if (land && iland_det > 0)  out[ichan][iloc] = 1.0;
           else if (snow && isnow_det > 0)  out[ichan][iloc] = 1.0;
           else if (ice && iice_det > 0)  out[ichan][iloc] = 1.0;
           else if (mixed && imix_det > 0)  out[ichan][iloc] = 1.0;
         }
       }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SurfTypeCheckRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
