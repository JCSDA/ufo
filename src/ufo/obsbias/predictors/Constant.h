/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_PREDICTORS_CONSTANT_H_
#define UFO_OBSBIAS_PREDICTORS_CONSTANT_H_

#include <string>
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

class Constant : public PredictorBase {
 public:
  explicit Constant(const eckit::Configuration &);
  ~Constant() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const std::vector<int> &,
               Eigen::MatrixXd &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_PREDICTORS_CONSTANT_H_
