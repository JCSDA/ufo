/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_COSINEOFLATITUDETIMESORBITNODE_H_
#define UFO_PREDICTORS_COSINEOFLATITUDETIMESORBITNODE_H_

#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class CosineOfLatitudeTimesOrbitNode : public PredictorBase {
 public:
  CosineOfLatitudeTimesOrbitNode(const eckit::Configuration &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  // default preconditioner for bias terms
  double precond_ = 0.01;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_COSINEOFLATITUDETIMESORBITNODE_H_
