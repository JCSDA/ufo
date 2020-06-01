/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ChannelUseflagCheckRad.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ChannelUseflagCheckRad>
       makerChannelUseflagCheckRad_("ChannelUseflagCheckRad");

// -----------------------------------------------------------------------------

ChannelUseflagCheckRad::ChannelUseflagCheckRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get channel use flags from options
  std::vector<int> useflag = options_.useflagChannel.value();
  ASSERT(useflag.size() == channels_.size());

  // TODO(EL) the following two lines will be removed when the revised filter behavior is in place
  // Include required variables from ObsDiag (note: included here to trigger posterFilter)
  invars_ += Variable("brightness_temperature_jacobian_surface_temperature@ObsDiag", channels_);
}

// -----------------------------------------------------------------------------

ChannelUseflagCheckRad::~ChannelUseflagCheckRad() {}

// -----------------------------------------------------------------------------

void ChannelUseflagCheckRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t nchans = channels_.size();

  // Get channel use flags from options
  std::vector<int> useflag = options_.useflagChannel.value();

  // Output
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ichan][iloc] = useflag[ichan];
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ChannelUseflagCheckRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
