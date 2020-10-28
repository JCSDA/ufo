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

Constant::Constant(const eckit::Configuration & conf, const std::vector<int> & jobs)
  : PredictorBase(conf, jobs) {
}

// -----------------------------------------------------------------------------

void Constant::compute(const ioda::ObsSpace & odb,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       ioda::ObsVector & out) const {
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  const std::size_t njobs = jobs_.size();
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    for (std::size_t jb = 0; jb < njobs; ++jb) {
      out[jl*njobs+jb] = 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
