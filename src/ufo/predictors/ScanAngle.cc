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
#include "ioda/ObsVector.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ScanAngle> makerFuncScanAngle_("sensorScanAngle");

// -----------------------------------------------------------------------------

ScanAngle::ScanAngle(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1)),
    var_name_(parameters.varName) {
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between scan_angle predictors of different orders
    name() = name() + "_order_" + std::to_string(order_);
  }
}

// -----------------------------------------------------------------------------

void ScanAngle::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        const ObsBias &,
                        ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  std::vector<float> view_angle(nlocs, 0.0);

  // retrieve the sensor view angle

  if (odb.dtype("MetaData", var_name_) == ioda::ObsDtype::Integer) {
    std::vector<int> view_angle2(nlocs, 0);
    odb.get_db("MetaData", var_name_, view_angle2);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      view_angle[jloc] = view_angle2[jloc]*1.0f;
    }
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
