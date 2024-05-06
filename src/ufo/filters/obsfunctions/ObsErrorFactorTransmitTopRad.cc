/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorTransmitTopRad.h"

#include <math.h>

#include <algorithm>
#include <limits>
#include <set>
#include <string>
#include <vector>
#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorTransmitTopRad>
       makerObsErrorFactorTransmitTopRad_("ObsErrorFactorTransmitTopRad");

// -----------------------------------------------------------------------------

ObsErrorFactorTransmitTopRad::ObsErrorFactorTransmitTopRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Include required variables from ObsDiag
  invars_ += Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_);
}

// -----------------------------------------------------------------------------

ObsErrorFactorTransmitTopRad::~ObsErrorFactorTransmitTopRad() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorTransmitTopRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ObsErrorFactorTransmitTopRad::compute start" << std::endl;
  // Get dimensions
  size_t nlocs = in.nlocs();
  if (nlocs == 0) {
    return;
  }
  size_t nchans = channels_.size();
  constexpr float eps = std::numeric_limits<float>::epsilon();
  // Inflate obs error as a function of model top-to-spaec transmittance
  std::vector<float> tao_top(nlocs);
  for (size_t ich = 0; ich < nchans; ++ich) {
    in.get(Variable("ObsDiag/transmittances_of_atmosphere_layer", channels_)[ich], 0, tao_top);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ich][iloc] = sqrt(1.0f / (std::abs(tao_top[iloc])+eps) );
    }
  }
  oops::Log::trace() << "ObsErrorFactorTransmitTopRad::compute done" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorTransmitTopRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
