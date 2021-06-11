/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_SCANANGLE_H_
#define UFO_PREDICTORS_SCANANGLE_H_
#include <string>
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

class ScanAngle : public PredictorBase {
 public:
  ScanAngle(const eckit::Configuration &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  std::string var_name_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_SCANANGLE_H_
