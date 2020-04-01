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

ObsFunctionObsErrorMean::ObsFunctionObsErrorMean(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Check required variables
  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.chList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get channels for CLW retrieval from options
  const std::vector<int> channels_clwret_ = {options_.ch238.value(), options_.ch314.value()};
  ASSERT(options_.ch238 !=0 && options_.ch314 !=0 && channels_clwret_.size() == 2);

  // Get variable groups for CLW retrieval
  ASSERT(options_.varGrp.value().size() == 2);

  // Get observation error model parameters from options
  ASSERT(channels_.size() == options_.clwClr.value().size());
  ASSERT(channels_.size() == options_.clwCld.value().size());
  ASSERT(channels_.size() == options_.obserrMin.value().size());
  ASSERT(channels_.size() == options_.obserrMax.value().size());

  // Get required variables from function
  conf_.set("clwret_ch238", options_.ch238);
  conf_.set("clwret_ch314", options_.ch314);
  conf_.set("clwret_types", options_.varGrp);
  conf_.set("test_group", options_.testGrp);
  conf_.set("bias_application", options_.addBias);
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

  // Get Mean CLW retrievals from function
  std::vector<float> clwretmean(nlocs);
  in.get(Variable("CLWRetMean@ObsFunction", conf_), clwretmean);

  const std::vector<float> &clw_clr_ = options_.clwClr;
  const std::vector<float> &clw_cld_ = options_.clwCld;
  const std::vector<float> &obserr_min_ = options_.obserrMin;
  const std::vector<float> &obserr_max_ = options_.obserrMax;
  // Calculate observation error for each channel
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      // constant obs error for clear sky
      if (clwretmean[iloc] <= clw_clr_[ichan]) {
        out[ichan][iloc] = obserr_min_[ichan];
      // observation error grows linearly with increaseing cloud amount under cloudy condition
      } else if (clwretmean[iloc] > clw_clr_[ichan] && clwretmean[iloc] < clw_cld_[ichan]) {
        float slope = (obserr_max_[ichan] - obserr_min_[ichan]) /
                      (clw_cld_[ichan] - clw_clr_[ichan]);
        out[ichan][iloc] = obserr_min_[ichan] + slope * (clwretmean[iloc] - clw_clr_[ichan]);
      // maxumum obs error for cloudy sky
      } else {
        out[ichan][iloc] = obserr_max_[ichan];
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
