/*
 * (C) Copyright 2020-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <cmath>
#include <string>
#include <vector>

#include "ufo/predictors/LapseRate.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<LapseRate> makerFuncLapseRate_("lapseRate");

// -----------------------------------------------------------------------------

LapseRate::LapseRate(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    lapseHelper_(parameters.lapseMean, vars),
    order_(parameters.order.value().value_or(1))
{
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between lapse_rate predictors of different orders
    name() = name() + "_order_" + std::to_string(order_);
  }

  if (vars.size() > 0) {
    geovars_ += LapseRateHelper::requiredGeovars();
    hdiags_ += lapseHelper_.requiredHdiags();
  } else {
    oops::Log::error() << "Channels size is ZERO !" << std::endl;
    ABORT("Channels size is ZERO !");
  }
}

// -----------------------------------------------------------------------------

void LapseRate::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs & geovals,
                        const ObsDiagnostics & ydiags,
                        const ObsBias &,
                        ioda::ObsVector & out) const {
  const std::size_t nvars = out.nvars();
  const std::size_t nlocs = out.nlocs();

  std::vector<std::vector<float>> tlapchn;
  std::vector<float> means;
  lapseHelper_.fillTlapchnAndMeans(odb, geovals, ydiags, tlapchn, means);

  if (nvars != tlapchn.size()) {
    oops::Log::error() << "LapseRate: nvars (" << nvars
                       << ") != number of channels (" << tlapchn.size() << ")"
                       << std::endl;
    ABORT("LapseRate: nvars mismatch");
  }

  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();
  const auto isBad = [fmiss](float v) { return v == fmiss || std::isnan(v); };

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      const float tlapchn_val = tlapchn[jvar][jloc];
      if (isBad(tlapchn_val) || isBad(means[jvar])) {
        out[jloc*nvars+jvar] = dmiss;
      } else {
        out[jloc*nvars+jvar] = std::pow((tlapchn_val - means[jvar]), order_);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
