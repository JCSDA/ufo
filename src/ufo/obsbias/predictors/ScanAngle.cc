/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/obsbias/predictors/ScanAngle.h"

#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ScanAngle> makerFuncScanAngle_("scan_angle");

// -----------------------------------------------------------------------------

ScanAngle::ScanAngle(const eckit::Configuration & conf)
  : PredictorBase(conf), order_(1) {
  // get the order if it is provided in options
  if (conf.has("predictor.options.order")) {
    conf.get("predictor.options.order", order_);

    // override the predictor name for differentiable
    name() = name() + "_order_" + std::to_string(order_);
  }
}

// -----------------------------------------------------------------------------

void ScanAngle::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        const std::vector<int> & jobs,
                        Eigen::MatrixXd & out) const {
  const std::size_t njobs = jobs.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.rows() == njobs && out.cols() == nlocs);

  // retrieve the sensor view angle
  std::vector<float> view_angle(nlocs, 0.0);
  odb.get_db("MetaData", "sensor_view_angle", view_angle);

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    out.col(jl).setConstant(pow(view_angle[jl]*Constants::deg2rad, order_));
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
