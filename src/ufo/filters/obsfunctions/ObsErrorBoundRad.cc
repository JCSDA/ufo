/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorBoundRad.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/obsfunctions/ObsErrorFactorLatRad.h"
#include "ufo/filters/obsfunctions/ObsErrorFactorTransmitTopRad.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorBoundRad> makerObsErrorBoundRad_("ObsErrorBoundRad");

// -----------------------------------------------------------------------------

ObsErrorBoundRad::ObsErrorBoundRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  const Variable &obserrlat = options_.obserrBoundLat.value();
  invars_ += obserrlat;

  const Variable &obserrtaotop = options_.obserrBoundTransmittop.value();

  invars_ += obserrtaotop;

  // Get test groups from options
  const std::string &errgrp_ = options_.testObserr.value();

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@"+errgrp_, channels_);
  invars_ += Variable("brightness_temperature@ObsError", channels_);
  invars_ += Variable("latitude@MetaData");
  invars_ += Variable("longitude@MetaData");
}

// -----------------------------------------------------------------------------

ObsErrorBoundRad::~ObsErrorBoundRad() {}

// -----------------------------------------------------------------------------

void ObsErrorBoundRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get observation error bounds from options
  const std::vector<float> &obserr_bound_max = options_.obserrBoundMax.value();
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get error factor from ObsFunction
  const Variable &obserrlat = options_.obserrBoundLat.value();
  ioda::ObsDataVector<float> errflat(in.obsspace(), obserrlat.toOopsVariables());
  in.get(obserrlat, errflat);

  // Get error factor from ObsFunction
  const Variable &obserrtaotop = options_.obserrBoundTransmittop.value();
  ioda::ObsDataVector<float> errftaotop(in.obsspace(), obserrtaotop.toOopsVariables());
  in.get(obserrtaotop, errftaotop);

  // Output integrated error bound for gross check
  std::vector<float> obserr(nlocs);
  std::vector<float> obserrdata(nlocs);
  const std::string &errgrp_ = options_.testObserr.value();
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+errgrp_, channels_)[ichan], obserrdata);
    in.get(Variable("brightness_temperature@ObsError", channels_)[ichan], obserr);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      float varinv = 1.0 / pow(obserrdata[iloc], 2);
      out[ichan][iloc] = obserr[iloc];
      if (varinv > 0.0) {
        out[ichan][iloc] = std::fmin(3.0 * obserr[iloc]
                               * (1.0 / pow(errflat[0][iloc], 2))
                               * (1.0 / pow(errftaotop[ichan][iloc], 2)), obserr_bound_max[ichan]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorBoundRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
