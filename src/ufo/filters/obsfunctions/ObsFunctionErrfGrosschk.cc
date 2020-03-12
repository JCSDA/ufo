/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionErrfGrosschk.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/obsfunctions/ObsFunctionErrfLat.h"
#include "ufo/filters/obsfunctions/ObsFunctionErrfTransmittop.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionErrfGrosschk> makerObsFuncErrfGrosschk_("ErrfGrosschk");

// -----------------------------------------------------------------------------

ObsFunctionErrfGrosschk::ObsFunctionErrfGrosschk(const eckit::LocalConfiguration conf)
  : invars_(), group_("ObsErrorData"), channels_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("channels") && conf_.has("obserr_max") && conf_.has("latitude_parameters"));

  // Check if using obserr from GSI for testing
  if (conf_.has("obserr_test")) group_ = conf_.getString("obserr_test");

  // Get channels from options
  const std::string chlist = conf.getString("channels");
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@"+group_, channels_);
  invars_ += Variable("brightness_temperature@ObsError", channels_);
  invars_ += Variable("latitude@MetaData");
  invars_ += Variable("longitude@MetaData");

  // Include required variables from ObsFunction
  ObsFunctionErrfTransmittop taotopfunc(conf_);
  invars_ += taotopfunc.requiredVariables();

  ObsFunctionErrfLat latfunc(conf_);
  invars_ += latfunc.requiredVariables();
}

// -----------------------------------------------------------------------------

ObsFunctionErrfGrosschk::~ObsFunctionErrfGrosschk() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrfGrosschk::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get observation error bounds from options
  std::vector<float> obserr_max = conf_.getFloatVector("obserr_max");

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get error factor from ObsFunction
  Variable fvar1(conf_);
  ioda::ObsDataVector<float> errftaotop(in.obsspace(), fvar1.toOopsVariables());
  ObsFunctionErrfTransmittop taotopfunc(conf_);
  taotopfunc.compute(in, errftaotop);

  // Get error factor from ObsFunction
  std::vector<float> errflat(nlocs);
  in.get(Variable("ErrfLat@ObsFunction", conf_), errflat);

  // Output integrated error bound for gross check
  std::vector<float> obserr(nlocs);
  std::vector<float> obserrdata(nlocs);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@"+group_, channels_)[ichan], obserrdata);
    in.get(Variable("brightness_temperature@ObsError", channels_)[ichan], obserr);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      float varinv = 1.0 / pow(obserrdata[iloc], 2);
      out[ichan][iloc] = obserr[iloc];
      if (varinv > 0.0) {
        out[ichan][iloc] = std::fmin(3.0f * obserr[iloc]
                               * (1.0f / pow(errflat[iloc], 2))
                               * (1.0f / pow(errftaotop[ichan][iloc], 2)), obserr_max[ichan]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionErrfGrosschk::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
