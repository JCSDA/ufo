/*
 * (C) Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/Logger.h"
#include "ufo/predictors/ReadBias.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ReadBias>
       makerFuncReadBias_("read_bias");

// -----------------------------------------------------------------------------

ReadBias::ReadBias(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
  group_name_(parameters.groupName) {
}

// -----------------------------------------------------------------------------

void ReadBias::compute(const ioda::ObsSpace & odb,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       const ObsBias &,
                       ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();
  const std::vector<std::string> variables = out.varnames().variables();

  out.zero();
  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    if (odb.has(group_name_, variables[jvar])) {
      std::vector<float> bias_value;
      odb.get_db(group_name_, variables[jvar], bias_value);
      for (std::size_t jloc = 0; jloc < nlocs; ++jloc)
        out[jloc*nvars+jvar] = bias_value[jloc];
    } else {
      const std::string errmsg = "Attempting to read bias from " + group_name_ + "/" +
                           variables[jvar] + " either the group or variable isn't in the ObsSpace.";
      throw eckit::BadParameter(errmsg, Here());
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
