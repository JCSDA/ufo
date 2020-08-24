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
                       ioda::ObsDataVector<double> & out) const {
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  for (auto & job : jobs_) {
    const std::string varname = name() + "_" + std::to_string(job);
    for (auto & item : out[varname]) {
      item = 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
