/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "ufo/predictors/Constant.h"

#include "ioda/ObsSpace.h"

namespace ufo {

static PredictorMaker<Constant> makerFuncConstant_("constant");

// -----------------------------------------------------------------------------

Constant::Constant(const eckit::Configuration & conf, const oops::Variables & vars)
  : PredictorBase(conf, vars) {
}

// -----------------------------------------------------------------------------

void Constant::compute(const ioda::ObsSpace & odb,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
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
