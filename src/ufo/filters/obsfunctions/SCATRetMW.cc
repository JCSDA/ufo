/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SCATRetMW.h"

#include <algorithm>
#include <cmath>
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

static ObsFunctionMaker<SCATRetMW> makerSCATRetMW_("SCATRetMW");

SCATRetMW::SCATRetMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Check required parameters
  // Get variable group types for scattering index retrieval from option
  ASSERT(options_.varGroup.value().size() == 1 || options_.varGroup.value().size() == 2);

  // Get channels for CLW retrieval from options
  const std::vector<int> channels_ = {options_.ch238.value(), options_.ch314.value(),
                                      options_.ch890.value()};
  ASSERT(options_.ch238 != 0 && options_.ch314 != 0 && options_.ch890 != 0 &&
         channels_.size() == 3);

  // Include list of required data from ObsSpace
  for (size_t igrp = 0; igrp < options_.varGroup.value().size(); ++igrp) {
    invars_ += Variable("brightness_temperature@" + options_.varGroup.value()[igrp], channels_);
  }
  invars_ += Variable("brightness_temperature@" + options_.testBias.value(), channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");
}

// -----------------------------------------------------------------------------

SCATRetMW::~SCATRetMW() {}

// -----------------------------------------------------------------------------

void SCATRetMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t ngrps = options_.varGroup.value().size();

  // Get required parameters
  const std::vector<std::string> &vargrp_ = options_.varGroup.value();
  const std::vector<int> channels_ = {options_.ch238.value(), options_.ch314.value(),
                                      options_.ch890.value()};

  // Get area fraction of each surface type from GeoVaLs
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Get observation, bias correction and HofX from ObsSpace
  std::vector<float> bt238(nlocs), bt314(nlocs), bt890(nlocs);
  for (size_t igrp = 0; igrp < ngrps; ++igrp) {
    // Get data based on group type
    in.get(Variable("brightness_temperature@"+vargrp_[igrp], channels_)[0], bt238);
    in.get(Variable("brightness_temperature@"+vargrp_[igrp], channels_)[1], bt314);
    in.get(Variable("brightness_temperature@"+vargrp_[igrp], channels_)[2], bt890);
    // Get bias based on group type
    if (options_.addBias.value() == vargrp_[igrp]) {
      std::vector<float> bias238(nlocs), bias314(nlocs), bias890(nlocs);
      in.get(Variable("brightness_temperature@"+options_.testBias.value(), channels_)[0], bias238);
      in.get(Variable("brightness_temperature@"+options_.testBias.value(), channels_)[1], bias314);
      in.get(Variable("brightness_temperature@"+options_.testBias.value(), channels_)[2], bias890);
      // Add bias correction to the assigned group
      if (options_.addBias.value() == "ObsValue") {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          bt238[iloc] = bt238[iloc] - bias238[iloc];
          bt314[iloc] = bt314[iloc] - bias314[iloc];
          bt890[iloc] = bt890[iloc] - bias890[iloc];
        }
      } else {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          bt238[iloc] = bt238[iloc] + bias238[iloc];
          bt314[iloc] = bt314[iloc] + bias314[iloc];
          bt890[iloc] = bt890[iloc] + bias890[iloc];
        }
      }
    }
    // Retrieve scattering index
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (water_frac[iloc] >= 0.99) {
        out[igrp][iloc] = -113.2 + (2.41 - 0.0049 * bt238[iloc]) * bt238[iloc]
                               + 0.454 * bt314[iloc] - bt890[iloc];
        out[igrp][iloc] = std::max(0.f, out[igrp][iloc]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SCATRetMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
