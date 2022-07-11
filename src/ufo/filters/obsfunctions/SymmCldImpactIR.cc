/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SymmCldImpactIR.h"

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

static ObsFunctionMaker<SymmCldImpactIR> makerSCIIR_("SymmCldImpactIR");

// -----------------------------------------------------------------------------

SymmCldImpactIR::SymmCldImpactIR(const eckit::LocalConfiguration config)
  : invars_(), channels_() {
  // Initialize options
  options_.deserialize(config);

  // Get channels from options
  std::string chlist = options_.chlist;
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include required variables from ObsDiag
  invars_ += Variable("brightness_temperature_assuming_clear_sky@ObsDiag", channels_);
}

// -----------------------------------------------------------------------------

SymmCldImpactIR::~SymmCldImpactIR() {}

// -----------------------------------------------------------------------------

void SymmCldImpactIR::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & SCI) const {
  const float missing = util::missingValue(missing);

  // Get dimensions
  size_t nlocs = in.nlocs();

  // Allocate vectors common across channels
  std::vector<float> clr(nlocs);
  std::vector<float> bak(nlocs);
  std::vector<float> obs(nlocs);
  std::vector<float> bias(nlocs);

  float Cmod, Cobs;

  for (size_t ich = 0; ich < SCI.nvars(); ++ich) {
    // Get channel-specific clr, bak, obs, and bias
    in.get(Variable("brightness_temperature_assuming_clear_sky@ObsDiag", channels_)[ich], clr);
    in.get(Variable("brightness_temperature@HofX", channels_)[ich], bak);
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ich], obs);
    if (in.has(Variable("brightness_temperature@ObsBiasData", channels_)[ich])) {
      in.get(Variable("brightness_temperature@ObsBiasData", channels_)[ich], bias);
    } else {
      std::fill(bias.begin(), bias.end(), 0.0f);
    }
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (clr[iloc] != missing && bak[iloc] != missing &&
          obs[iloc] != missing && bias[iloc] != missing) {
        // Temporarily account for ZERO clear-sky BT output from CRTM
        // TODO(JJG): change CRTM clear-sky behavior
        if (clr[iloc] > -1.0f && clr[iloc] < 1.0f) clr[iloc] = bak[iloc];

        // HofX contains bias correction; subtracting it here
        Cmod = std::abs(clr[iloc] - bak[iloc] + bias[iloc]);
        Cobs = std::abs(clr[iloc] - obs[iloc] + bias[iloc]);
        SCI[ich][iloc] = 0.5f * (Cmod + Cobs);
      } else {
        SCI[ich][iloc] = missing;
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SymmCldImpactIR::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
