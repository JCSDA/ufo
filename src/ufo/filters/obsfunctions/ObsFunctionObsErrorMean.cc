/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionObsErrorMean.h"

#include <math.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/obsfunctions/ObsFunctionCLWRetMean.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionObsErrorMean> makerObsFuncObsErrorMean_("ObsErrorMean");

// ----------------------------------------------------------------------------------

ObsFunctionObsErrorMean::ObsFunctionObsErrorMean(const eckit::LocalConfiguration conf)
  : invars_(), channels_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("clwret_type") && conf_.has("channels") &&
         conf_.has("clw_clr") && conf_.has("clw_cld") &&
         conf_.has("obserr_clr") && conf_.has("obserr_cld"));

  // Get channels from options
  const std::string chlist = conf_.getString("channels");
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Get required variables from function
  ObsFunctionCLWRetMean clwretfunc(conf_);
  invars_ += clwretfunc.requiredVariables();
}

// -----------------------------------------------------------------------------

ObsFunctionObsErrorMean::~ObsFunctionObsErrorMean() {}

// -----------------------------------------------------------------------------

void ObsFunctionObsErrorMean::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get parameters for observation errors from options
  std::vector<float> clw_clr = conf_.getFloatVector("clw_clr");
  std::vector<float> clw_cld = conf_.getFloatVector("clw_cld");
  std::vector<float> obserr_clr = conf_.getFloatVector("obserr_clr");
  std::vector<float> obserr_cld = conf_.getFloatVector("obserr_cld");

  // Get Mean CLW retrievals from function
  std::vector<float> clwretmean(nlocs);
  in.get(Variable("CLWRetMean@ObsFunction", conf_), clwretmean);

  // Calculate observation error for each channel
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (clwretmean[iloc] <= clw_clr[ichan]) {
        out[ichan][iloc] = obserr_clr[ichan];
      } else if (clwretmean[iloc] > clw_clr[ichan] && clwretmean[iloc] < clw_cld[ichan]) {
        out[ichan][iloc] = obserr_clr[ichan] + (clwretmean[iloc] - clw_clr[ichan]) *
                          (obserr_cld[ichan] - obserr_clr[ichan]) /
                          (clw_cld[ichan] - clw_clr[ichan]);
      } else {
        out[ichan][iloc] = obserr_cld[ichan];
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionObsErrorMean::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
