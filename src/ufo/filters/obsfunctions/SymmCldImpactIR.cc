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
#include <sstream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<SymmCldImpactIR> makerSCIIR_("SymmCldImpactIR");

// -----------------------------------------------------------------------------

SymmCldImpactIR::SymmCldImpactIR(const eckit::LocalConfiguration config)
  : invars_(), channels_() {
  oops::Log::debug() << "SymmCldImpactIR: config = " << config << std::endl;
  // Initialize options
  options_.deserialize(config);

  // Get channels from options
  std::string chlist = options_.chlist;
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include required variables from ObsDiag
  invars_ += Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_);
}

// -----------------------------------------------------------------------------

SymmCldImpactIR::~SymmCldImpactIR() {}

// -----------------------------------------------------------------------------

void SymmCldImpactIR::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & SCI) const {
  const float missing = util::missingValue<float>();

  // Get dimensions
  size_t nlocs = in.nlocs();

  // Is the option for scaling by obs-minus-background enabled
  bool scale_by_omb = options_.scale_by_omb.value();

  // Allocate vectors common across channels
  std::vector<float> clr(nlocs);
  std::vector<float> bak(nlocs);
  std::vector<float> obs(nlocs);
  std::vector<float> bias(nlocs);

  const float c1 = options_.sigmoid_c1.value();    // Affects slope of sigmoid function
  const float c2 = options_.sigmoid_c2.value();    // The 50-percent value of sigmoid function

  float Cmod, Cobs, Comb, dx, frac;

  for (size_t ich = 0; ich < SCI.nvars(); ++ich) {
    // Get channel-specific clr, bak, obs, and bias
    in.get(Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_)[ich], clr);
    in.get(Variable("HofX/brightnessTemperature", channels_)[ich], bak);
    in.get(Variable("ObsValue/brightnessTemperature", channels_)[ich], obs);
    if (in.has(Variable("ObsBiasData/brightnessTemperature", channels_)[ich])) {
      in.get(Variable("ObsBiasData/brightnessTemperature", channels_)[ich], bias);
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
        if (scale_by_omb) {
          Comb = std::min(std::abs(obs[iloc] - bak[iloc]), 100.0f);
          dx = (Comb - c2)*0.01f;
          frac = std::max(0.1f, 1.0f/(1.0f + exp(-c1*dx)));
          SCI[ich][iloc] = frac * SCI[ich][iloc];
        }
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
