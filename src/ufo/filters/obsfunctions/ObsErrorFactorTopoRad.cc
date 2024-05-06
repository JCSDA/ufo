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
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorTopoRad>
       makerObsErrorFactorTopoRad_("ObsErrorFactorTopoRad");

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::ObsErrorFactorTopoRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);
  ASSERT(inst == "amsua" || inst == "atms"    || inst == "mhs" ||
         inst == "iasi" || inst == "cris-fsr" || inst == "airs" ||
         inst == "avhrr3" || inst == "seviri");

  if (inst == "amsua" || inst == "atms") {
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
  invars_ += Variable("GeoVaLs/surface_geopotential_height");
}

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::~ObsErrorFactorTopoRad() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorTopoRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_)[0]);

  // Get surface geopotential height
  std::vector<float> zsges(nlocs);
  in.get(Variable("GeoVaLs/surface_geopotential_height"), zsges);

  // Inflate obs error as a function of terrian height (>2000) and surface-to-space transmittance
  if (inst == "iasi" || inst == "cris-fsr" || inst == "airs" || inst == "avhrr3" ||
      inst == "seviri") {
    std::vector<float> tao_sfc(nlocs);
    for (size_t ich = 0; ich < nchans; ++ich) {
      in.get(Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_)[ich],
             nlevs - 1, tao_sfc);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ich][iloc] = 1.0;
        if (zsges[iloc] > 2000.0) {
          float factor = pow((2000.0/zsges[iloc]), 4);
          out[ich][iloc] = sqrt(1.0 / (1.0 - (1.0 - factor) * tao_sfc[iloc]));
        }
      }
    }
  } else if (inst == "mhs") {
    std::vector<float> tao_sfc(nlocs);
    for (size_t ich = 0; ich < nchans; ++ich) {
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ich][iloc] = 1.0;
        if (zsges[iloc] > 2000.0) {
          float factor = 2000.0/zsges[iloc];
          out[ich][iloc] = sqrt(1.0 / factor);
        }
      }
    }
  } else if (inst == "amsua" || inst == "atms") {
    // Set channel numbers
    int ich238, ich314, ich503, ich528, ich536, ich544, ich549, ich890;
    if (inst == "amsua") {
      ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 4, ich536 = 5;
      ich544 = 6, ich549 = 7, ich890 = 15;
    } else if (inst == "atms") {
      ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 5, ich536 = 6;
      ich544 = 7, ich549 = 8, ich890 = 16;
    }

    float factor;
    std::vector<int> qcflagdata;
    std::vector<float> obserrdata;
    const std::string &errgrp = options_.testObserr.value();
    const std::string &flaggrp = options_.testQCflag.value();
    const float missing = util::missingValue<float>();

    // Calculate error factors (error_factors) for each channel
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      size_t channel = ichan + 1;
      in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
      in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ichan][iloc] = 1.0;
        if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                             : qcflagdata[iloc] = 0;
        (qcflagdata[iloc] != 0) ? (factor = 0.0) : (factor = 1.0);

        if (zsges[iloc] > 2000.0) {
          if (channel <= ich544 || channel >= ich890) {
            out[ichan][iloc] = (2000.0/zsges[iloc]) * factor;
          }
          if ((zsges[iloc] > 4000.0) && (channel == ich549)) {
            out[ichan][iloc] = (4000.0/zsges[iloc]) * factor;
          }
          if (factor > 0.0) out[ichan][iloc] = sqrt(1.0 / out[ichan][iloc]);
        }
      }
    }
  } else {
    oops::Log::error() << "ObsErrorFactorTopoRad: Invalid instrument (sensor) specified: " << inst
                       << "  The valid instruments are: iasi, cris-fsr, airs, avhrr3, seviri, "
                       << "  amsua, atms, and mhs."
                       << std::endl;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorTopoRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
