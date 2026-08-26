/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorTopoRad.h"

#include <math.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorTopoRad>
       makerObsErrorFactorTopoRad_("ObsErrorFactorTopoRad");

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::ObsErrorFactorTopoRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorTopoRad constructor" << std::endl;
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Build the channel -> height map.
  // If no height groups are specified, apply the default reference height (2000 m) to every
  // channel. Only channels listed in a height group are corrected, using the listed height.
  const auto &heightGroups = options_.heightGroups.value();
  if (heightGroups.empty()) {
    for (int channel : channels_) {
      channelHeights_[channel] = 2000.0;
    }
  } else {
    for (const auto &group : heightGroups) {
      std::set<int> groupChannels = oops::parseIntSet(group.channelList);
      for (int channel : groupChannels) {
        channelHeights_[channel] = group.height;
      }
    }
  }

  // Get sensor information from options
  sensor_ = options_.sensor.value();
  ASSERT(sensor_ == "infrared" || sensor_ == "microwave");

  if (sensor_ == "infrared") {
    // Get test groups from options
    const std::string &errgrp = options_.testObserr.value();
    const std::string &flaggrp = options_.testQCflag.value();

    // Include list of required data from ObsSpace
    invars_ += Variable(errgrp+"/brightnessTemperature", channels_);
    invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);
  }

  // Include required variables from ObsDiag
  invars_ += Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/geopotential_height_at_surface");
}

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::~ObsErrorFactorTopoRad() {
  oops::Log::trace() << "ObsErrorFactorTopoRad destructor"  << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorTopoRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ObsErrorFactorTopoRad compute start" << std::endl;

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_)[0]);

  // Get surface geopotential height
  std::vector<float> zsges(nlocs);
  in.get(Variable("GeoVaLs/geopotential_height_at_surface"), zsges);

  // Inflate obs error as a function of terrain height and surface-to-space transmittance
  if (sensor_ == "infrared") {
    std::vector<float> tao_sfc(nlocs);
    for (size_t ich = 0; ich < nchans; ++ich) {
      in.get(Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_)[ich],
             nlevs - 1, tao_sfc);
      const auto it = channelHeights_.find(channels_[ich]);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ich][iloc] = 1.0;
        if (it != channelHeights_.end() && zsges[iloc] > it->second) {
          float factor = pow((it->second/zsges[iloc]), 4);
          out[ich][iloc] = sqrt(1.0 / (1.0 - (1.0 - factor) * tao_sfc[iloc]));
        }
      }
    }
  } else if (sensor_ == "microwave") {
    float factor;
    std::vector<int> qcflagdata;
    std::vector<float> obserrdata;
    const std::string &errgrp = options_.testObserr.value();
    const std::string &flaggrp = options_.testQCflag.value();
    const float missing = util::missingValue<float>();

    for (size_t ich = 0; ich < nchans; ++ich) {
      in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ich], obserrdata);
      in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ich], qcflagdata);
      const auto it = channelHeights_.find(channels_[ich]);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ich][iloc] = 1.0;
        if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                             : qcflagdata[iloc] = 0;
        (qcflagdata[iloc] != 0) ? (factor = 0.0) : (factor = 1.0);

        if (it != channelHeights_.end() && zsges[iloc] > it->second) {
          out[ich][iloc] = it->second/zsges[iloc] * factor;
        }
        if (factor > 0.0) {
          out[ich][iloc] = sqrt(1.0 / out[ich][iloc]);
        }
      }
    }
  } else {  // Already checked but just to include an else statement here
    oops::Log::error() << "ObsErrorFactorTopoRad: Invalid sensor specified: " << sensor_
                       << "  The valid sensors are: infrared and microwave."
                       << std::endl;
  }
  oops::Log::trace() << "ObsErrorFactorTopoRad compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorTopoRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
