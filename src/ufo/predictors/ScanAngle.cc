/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/predictors/ScanAngle.h"

#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ScanAngle> makerFuncScanAngle_("scan_angle");

// -----------------------------------------------------------------------------

ScanAngle::ScanAngle(const eckit::Configuration & conf, const oops::Variables & vars)
  : PredictorBase(conf, vars), order_(1) {
  // get the order if it is provided in options
  if (conf.has("options.order")) {
    conf.get("options.order", order_);

    // override the predictor name for differentiable
    name() = name() + "_order_" + std::to_string(order_);
  }

  if (conf.has("options.var_name")) {
    conf.get("options.var_name", var_name_);
  }
}

// -----------------------------------------------------------------------------

void ScanAngle::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  // retrieve the sensor view angle
  std::vector<float> view_angle(nlocs, 0.0);
  if ( var_name_.empty() ) {
    odb.get_db("MetaData", "sensor_view_angle", view_angle);
  } else {
    odb.get_db("MetaData", var_name_, view_angle);
  }

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = pow(view_angle[jloc] * Constants::deg2rad, order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
