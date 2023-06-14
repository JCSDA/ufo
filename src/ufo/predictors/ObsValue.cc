/*
 * (C) Copyright 2023 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/predictors/ObsValue.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ObsValue>
       makerFuncObsValue_("obsvalue");

// -----------------------------------------------------------------------------

ObsValue::ObsValue(const Parameters_ & parameters, const oops::Variables & vars)
  : PredictorBase(parameters, vars),
  var_name_(parameters.varName) {
}


// -----------------------------------------------------------------------------

void ObsValue::compute(const ioda::ObsSpace & odb,
                             const GeoVaLs &,
                             const ObsDiagnostics &,
                             const ObsBias &,
                             ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();

  std::vector<float> obs_value(nlocs, 0.0);
  odb.get_db("ObsValue", var_name_, obs_value);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = obs_value[jloc];
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
