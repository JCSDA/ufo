/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorBoundIR.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorBoundIR> makerObsErrorBoundIR_("ObsErrorBoundIR");

// -----------------------------------------------------------------------------

ObsErrorBoundIR::ObsErrorBoundIR(const eckit::LocalConfiguration & conf)
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
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();

  // Include list of required data from ObsSpace
  invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);
  invars_ += Variable(errgrp+"/brightnessTemperature", channels_);
}

// -----------------------------------------------------------------------------

ObsErrorBoundIR::~ObsErrorBoundIR() {}

// -----------------------------------------------------------------------------

void ObsErrorBoundIR::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ObsErrorBoundIR::compute start" << std::endl;
  // Get observation error bounds from options
  const std::vector<float> &obserr_bound_max = options_.obserrBoundMax.value();
  // Get dimensions
  size_t nlocs = in.nlocs();
  if (nlocs == 0) {
    return;
  }
  size_t nchans = channels_.size();

  // Get error factor from ObsFunction
  const Variable &obserrlat = options_.obserrBoundLat.value();
  ioda::ObsDataVector<float> errflat(in.obsspace(), obserrlat.toOopsObsVariables());
  in.get(obserrlat, errflat);

  // Get error factor from ObsFunction
  const Variable &obserrtaotop = options_.obserrBoundTransmittop.value();
  ioda::ObsDataVector<float> errftaotop(in.obsspace(), obserrtaotop.toOopsObsVariables());
  in.get(obserrtaotop, errftaotop);

  // Get original observation error.  If not explicitly passed through the YAML, check the ObsSpace
  // This channel-dependent error is assumed to be constant across all obs locations.
  std::vector<float> obserr(nchans, 0.0f);
  if (options_.obserrOriginal.value() != boost::none) {
    obserr = options_.obserrOriginal.value().get();
  } else {
  // Get original observation error (uninflated) from ObsSpace
      std::vector<std::vector<float>> obserr2(nchans, std::vector<float>(nlocs));
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      in.get(Variable("ObsError/brightnessTemperature", channels_)[ichan], obserr2[ichan]);
      obserr[ichan] = obserr2[ichan][0];
    }
  }

  // Output integrated error bound for gross check
  std::vector<float> obserrdata(nlocs);  //!< effective obs err
  std::vector<int> qcflagdata(nlocs);    //!< effective qcflag
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const float missing = util::missingValue<float>();
  float varinv = 0.0;
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);
    in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                           : qcflagdata[iloc] = 0;
      (qcflagdata[iloc] == 0) ? (varinv = 1.0 / pow(obserrdata[iloc], 2)) : (varinv = 0.0);
      out[ichan][iloc] = obserr[ichan];
      if (varinv > 0.0) {
        out[ichan][iloc] = std::fmin(3.0 * obserr[ichan]
                               * (1.0 / pow(errflat[0][iloc], 2))
                               * (1.0 / pow(errftaotop[ichan][iloc], 2)), obserr_bound_max[ichan]);
      }
    }
  }
  oops::Log::trace() << "ObsErrorBoundIR::compute end" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorBoundIR::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
