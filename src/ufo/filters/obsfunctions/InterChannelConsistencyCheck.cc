/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/InterChannelConsistencyCheck.h"

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
#include "ufo/utils/StringUtils.h"

namespace ufo {

static ObsFunctionMaker<InterChannelConsistencyCheck>
       makerInterChannelConsistencyCheck_("InterChannelConsistencyCheck");

// -----------------------------------------------------------------------------

InterChannelConsistencyCheck::InterChannelConsistencyCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // TODO(EL) the following two lines will be removed when the revised filter behavior is in place
  // Include required variables from ObsDiag (note: included here to trigger posterFilter)
  invars_ += Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag", channels_);
}

// -----------------------------------------------------------------------------

InterChannelConsistencyCheck::~InterChannelConsistencyCheck() {}

// -----------------------------------------------------------------------------

void InterChannelConsistencyCheck::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);

  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t nchans = channels_.size();

  // Get channel use flags from options
  const std::vector<int> &use_flag = options_.useflagChannel.value();

  // Get test groups from options
  const std::string &flaggrp = options_.testQCflag.value();
  const std::string &errgrp = options_.testObserr.value();

  // Get effective observation error and qcflag from ObsSpace
  // Convert effective observation error to inverse of the error variance
  const float missing = util::missingValue(missing);
  std::vector<int> qcflagdata(nlocs);
  std::vector<float> obserrdata(nlocs);
  std::vector<std::vector<float>> varinv(nchans, std::vector<float>(nlocs, 0.0));
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+errgrp, channels_)[ichan], obserrdata);
    in.get(Variable("brightness_temperature@"+flaggrp, channels_)[ichan], qcflagdata);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                           : qcflagdata[iloc] = 0;
      (qcflagdata[iloc] == 0) ? (varinv[ichan][iloc] = 1.0 / pow(obserrdata[iloc], 2))
                              : (varinv[ichan][iloc] = 0.0);
    }
  }

  // Inter-channel consistency check
  const size_t ncheck = 6;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    for (int ichan = 0; ichan < nchans; ++ichan) out[ichan][iloc] = 0;
    int kval = 0;
    for (int ichan = 1; ichan < ncheck; ++ichan) {
      int channel = ichan + 1;
      if (varinv[ichan][iloc] <= 0.0 && use_flag[ichan] >= 1) {
        kval = std::max(channel-1, kval);
        if ((inst == "amsua" || inst == "atms") && channel <= 3) kval = 0;
      }
    }
    if (kval > 0) {
      for (int ichan = 0; ichan < kval; ++ichan) out[ichan][iloc] = 1;
      int channel = 15;
      out[channel-1][iloc] = 1;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & InterChannelConsistencyCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
