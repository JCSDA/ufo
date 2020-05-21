/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_PREDICTORS_SCANANGLE_H_
#define UFO_OBSBIAS_PREDICTORS_SCANANGLE_H_

#include <vector>

#include "ufo/obsbias/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class ScanAngle : public PredictorBase {
 public:
  ScanAngle(const eckit::Configuration &, const std::vector<int> &);
  ~ScanAngle() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               Eigen::MatrixXd &) const override;

 private:
  int order_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_PREDICTORS_SCANANGLE_H_
