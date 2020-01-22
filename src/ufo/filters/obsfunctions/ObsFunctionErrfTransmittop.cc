/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionErrfTransmittop.h"

#include <math.h>

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

static ObsFunctionMaker<ObsFunctionErrfTransmittop> makerObsFuncErrfTransmittop_("ErrfTransmittop");

// -----------------------------------------------------------------------------

ObsFunctionErrfTransmittop::ObsFunctionErrfTransmittop(const eckit::LocalConfiguration conf)
  : invars_(), channels_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("channels"));

  // Get channels from options
  const std::string chlist = conf_.getString("channels");
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include required variables from ObsDiag
  invars_ += Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_);
}

// -----------------------------------------------------------------------------

ObsFunctionErrfTransmittop::~ObsFunctionErrfTransmittop() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrfTransmittop::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Inflate obs error as a function of model top-to-spaec transmittance
  std::vector<float> tao_top(nlocs);
  for (size_t ich = 0; ich < nchans; ++ich) {
    in.get(Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_)[ich], 1, tao_top);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ich][iloc] = sqrt(1.0 / tao_top[iloc]);
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionErrfTransmittop::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
