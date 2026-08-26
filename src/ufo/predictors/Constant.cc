/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsVector.h"
#include "oops/util/abor1_cpp.h"
#include "ufo/GeoVaLs.h"
#include "ufo/predictors/Constant.h"

namespace ufo {

static PredictorMaker<Constant> makerFuncConstant_("constant");

// -----------------------------------------------------------------------------

Constant::Constant(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
  surface_(parameters.surface) {
  if (surface_ == "all") {
    name() = "constant";
  } else {
    geovars_ += oops::Variables({"landmask"});
    if (surface_ == "land only") {
      name() = name() + "_land";
    } else if (surface_ == "sea only") {
      name() = name() + "_sea";
    } else if (surface_ == "land sea mask") {
      name() = name() + "_slmask";
    }
  }
}

// -----------------------------------------------------------------------------

void Constant::compute(const ioda::ObsSpace &,
                       const GeoVaLs & geovals,
                       const ObsDiagnostics &,
                       const ObsBias &,
                       ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();

  out.zero();
  if (surface_ == "all") {
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        out[jloc*nvars+jvar] = 1.0;
      }
    }
  } else if (surface_ == "land only" || surface_ == "sea only" || surface_ == "land sea mask") {
    // get landmask only when needed
    std::vector<int> slmask(nlocs, 0.0);
    geovals.get(slmask, oops::Variable{"landmask"});

    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        out[jloc*nvars+jvar] = 1.0;
        // assign 0 to constant predictor over land
        if (surface_ == "sea only" && slmask[jloc] == 1) {
          out[jloc*nvars+jvar] = 0.0;
        }
        // assign 0 to constant predictor over sea
        if (surface_ == "land only" && slmask[jloc] == 0) {
          out[jloc*nvars+jvar] = 0.0;
        }
        // assign 1 to constant predictor over land and -1 over sea
        if (surface_ == "land sea mask" && slmask[jloc] == 0) {
            out[jloc*nvars+jvar] = -1.0;
        }
      }
    }
  } else {
    ABORT("Surface selection is not applicable !");
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
