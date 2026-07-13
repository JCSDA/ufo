/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SymmCldImpact.h"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/SymmCldImpactUtils.h"

namespace ufo {

static ObsFunctionMaker<SymmCldImpact> makerSCI_("SymmCldImpact");

// -----------------------------------------------------------------------------

SymmCldImpact::SymmCldImpact(const eckit::LocalConfiguration config)
  : invars_(), channels_() {
  oops::Log::trace() << "SymmCldImpact constructor" << std::endl;
  oops::Log::debug() << "SymmCldImpact: config = " << config << std::endl;

  options_.deserialize(config);

  std::string chlist = options_.chlist;
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  if (options_.btlim.value().has_value()) {
    const std::vector<float> & btlim = options_.btlim.value().value();
    if (btlim.size() != channels_.size()) {
      oops::Log::error() << "SymmCldImpact: btlim size (" << btlim.size()
                         << ") does not match number of channels ("
                         << channels_.size() << ")!" << std::endl;
      ABORT("SymmCldImpact: btlim size mismatch!");
    }
    // Harnish2016: only needs HofX brightness temperature
    invars_ += Variable("HofX/brightnessTemperature", channels_);
  } else {
    // Okamoto2014: needs clear-sky BT from ObsDiag and HofX
    invars_ += Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_);
    invars_ += Variable("HofX/brightnessTemperature", channels_);
  }
}

// -----------------------------------------------------------------------------

SymmCldImpact::~SymmCldImpact() {
  oops::Log::trace() << "SymmCldImpact destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void SymmCldImpact::compute(const ObsFilterData & in,
                             ioda::ObsDataVector<float> & SCI) const {
  oops::Log::trace() << "SymmCldImpact compute start" << std::endl;

  const float fmiss = util::missingValue<float>();
  const size_t nlocs = in.nlocs();
  const size_t nvars = SCI.nvars();

  // order_ not needed for obsfunction
  const int order = 1;

  std::vector<float> bak(nlocs);
  std::vector<float> obs(nlocs);
  std::vector<float> bias(nlocs, 0.0f);
  std::vector<float> out(nlocs * nvars, fmiss);

  if (options_.btlim.value().has_value()) {
    // Harnish2016
    const std::vector<float> btlim = options_.btlim.value().value();
    for (size_t ich = 0; ich < nvars; ++ich) {
      in.get(Variable("HofX/brightnessTemperature", channels_)[ich], bak);
      in.get(Variable("ObsValue/brightnessTemperature", channels_)[ich], obs);
      if (in.has(Variable("ObsBiasData/brightnessTemperature", channels_)[ich])) {
        in.get(Variable("ObsBiasData/brightnessTemperature", channels_)[ich], bias);
        // HofX contains the bias correction; subtract it here
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          if (bias[iloc] != fmiss) {
            bak[iloc] -= bias[iloc];
          } else {
            bak[iloc] = fmiss;
          }
        }
      }
      computeSCIHarnish(bak, obs, btlim[ich],
                        nlocs, ich, nvars, order, out);
    }
  } else {
    // Okamoto2014
    std::vector<float> clr(nlocs);
    const bool scale_by_omb = options_.scale_by_omb.value();
    const float c1 = options_.sigmoid_c1.value();
    const float c2 = options_.sigmoid_c2.value();
    for (size_t ich = 0; ich < nvars; ++ich) {
      in.get(Variable("ObsDiag/brightness_temperature_assuming_clear_sky",
                       channels_)[ich], clr);
      in.get(Variable("HofX/brightnessTemperature", channels_)[ich], bak);
      in.get(Variable("ObsValue/brightnessTemperature", channels_)[ich], obs);
      if (in.has(Variable("ObsBiasData/brightnessTemperature", channels_)[ich])) {
        in.get(Variable("ObsBiasData/brightnessTemperature", channels_)[ich], bias);
        // HofX contains the bias correction; subtract it here
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          if (bias[iloc] != fmiss) {
            bak[iloc] -= bias[iloc];
          } else {
            bak[iloc] = fmiss;
          }
        }
      }
      computeSCIOkamoto(clr, bak, obs, scale_by_omb, c1, c2,
                        nlocs, ich, nvars, order, out);
    }
  }

  // Copy output into SCI ObsDataVector
  for (size_t ich = 0; ich < nvars; ++ich) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      SCI[ich][iloc] = out[iloc * nvars + ich];
    }
  }

  oops::Log::trace() << "SymmCldImpact compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & SymmCldImpact::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
