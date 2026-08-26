/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/SymmCldImpactPred.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/SymmCldImpactUtils.h"

namespace ufo {

static PredictorMaker<SymmCldImpactPred>
       makerFuncSymmCldImpactPred_("SymmCldImpact");

// -----------------------------------------------------------------------------

SymmCldImpactPred::SymmCldImpactPred(const Parameters_ & parameters,
                                     const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1)) {
  oops::Log::trace() << "SymmCldImpactPred constructor" << std::endl;

  options_ = parameters;

  if (parameters.order.value() != boost::none) {
    // Override predictor name to distinguish between orders
    name() = name() + "_order_" + std::to_string(order_);
  }

  if (vars.size() > 0) {
    // Parse and validate channels for both formulations
    std::string chlist = options_.chlist;
    std::set<int> channelset = oops::parseIntSet(chlist);
    std::copy(channelset.begin(), channelset.end(),
              std::back_inserter(channels_));
    std::vector<int> varschans(vars.channels().begin(), vars.channels().end());
    std::sort(varschans.begin(), varschans.end());
    if (channels_ != varschans) {
      oops::Log::error() << "SymmCldImpactPred: channels list does not "
                         << "match obs operator channels!" << std::endl;
      ABORT("SymmCldImpactPred: channels mismatch!");
    }
    if (options_.btlim.value().has_value()) {
      // Harnish2016: needs HofX brightness temperature
      // In unit tests, use_hofx=true to avoid ObsDiagnostics initialization
      // crash in test_Predictor.x; in real DA runs, brightness_temperature
      // is provided by the forward operator via ObsDiagnostics
      const std::vector<float> & btlim = options_.btlim.value().value();
      if (btlim.size() != channels_.size()) {
        oops::Log::error() << "SymmCldImpactPred: btlim size (" << btlim.size()
                           << ") does not match number of channels ("
                           << channels_.size() << ")!" << std::endl;
        ABORT("SymmCldImpactPred: btlim size mismatch!");
      }
      if (!options_.use_hofx.value()) {
        hdiags_ += oops::ObsVariables({"brightness_temperature"}, vars.channels());
      }
    } else {
      // Okamoto2014: needs clear-sky and all-sky brightness temperature
      hdiags_ += oops::ObsVariables({"brightness_temperature_assuming_clear_sky"},
                                     vars.channels());
      hdiags_ += oops::ObsVariables({"brightness_temperature"}, vars.channels());
    }
  } else {
    oops::Log::error() << "SymmCldImpactPred: channels size is ZERO!" << std::endl;
    ABORT("SymmCldImpactPred: channels size is ZERO!");
  }
}
// -----------------------------------------------------------------------------

void SymmCldImpactPred::compute(const ioda::ObsSpace & odb,
                                const GeoVaLs &,
                                const ObsDiagnostics & ydiags,
                                const ObsBias &,
                                ioda::ObsVector & out) const {
  oops::Log::trace() << "SymmCldImpactPred compute start" << std::endl;

  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();
  const float fmiss = util::missingValue<float>();

  std::vector<float> bak(nlocs, 0.0f);
  std::vector<float> obs(nlocs, 0.0f);
  std::vector<float> outd(nlocs * nvars, fmiss);

  if (options_.btlim.value().has_value()) {
    // Harnish2016
    const std::vector<float> btlim = options_.btlim.value().value();
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      const std::string chstr = std::to_string(channels_[jvar]);
      // In unit tests (use_hofx=true), read from HofX in obs database;
      // in real DA runs, read from ObsDiagnostics provided by forward operator
      if (options_.use_hofx.value()) {
        odb.get_db("HofX", "brightnessTemperature", bak, {channels_[jvar]});
      } else {
        ydiags.get(bak, "brightness_temperature_" + chstr);
      }
      odb.get_db("ObsValue", "brightnessTemperature", obs, {channels_[jvar]});
      computeSCIHarnish(bak, obs, btlim[jvar],
                        nlocs, jvar, nvars, order_, outd);
    }
  } else {
    // Okamoto2014
    std::vector<float> clr(nlocs, 0.0f);
    const bool scale_by_omb = options_.scale_by_omb.value();
    const float c1 = options_.sigmoid_c1.value();
    const float c2 = options_.sigmoid_c2.value();
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      const std::string chstr = std::to_string(channels_[jvar]);
      ydiags.get(clr, "brightness_temperature_assuming_clear_sky_" + chstr);
      ydiags.get(bak, "brightness_temperature_" + chstr);
      odb.get_db("ObsValue", "brightnessTemperature", obs, {channels_[jvar]});
      computeSCIOkamoto(clr, bak, obs, scale_by_omb, c1, c2,
                        nlocs, jvar, nvars, order_, outd);
    }
  }

  // Copy to the output ObsVector, replacing the float missing value with the
  // double missing value expected by ObsBiasOperator.
  const double dmiss = util::missingValue<double>();
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      const float v = outd[jloc * nvars + jvar];
      out[jloc * nvars + jvar] = (v == fmiss) ? dmiss : static_cast<double>(v);
    }
  }

  oops::Log::trace() << "SymmCldImpactPred compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
