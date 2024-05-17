/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/Constant.h"

#include "ioda/ObsVector.h"

namespace ufo {

static PredictorMaker<Constant> makerFuncConstant_("constant");

// -----------------------------------------------------------------------------

Constant::Constant(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars) {
}

// -----------------------------------------------------------------------------

void Constant::compute(const ioda::ObsSpace &,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       const ObsBias &,
                       ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
