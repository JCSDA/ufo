/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_ORBITALANGLE_H_
#define UFO_PREDICTORS_ORBITALANGLE_H_
#include <string>
#include <vector>
#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------
/**
 *This orbital angle predictor is used to fit residual errors as a function of satellite
 *orbital angle using a Fourier series. The data must contain this orbital angle time
 *series in the variable "satellite_orbital_angle@MetaData". Two member variables are 
 *used to store the order of the term in the series being calculated (order_) and the 
 *Fourier component (cos or sin), respectively. These are read from the yaml configuration
 *file.
 */

class OrbitalAngle : public PredictorBase {
 public:
  OrbitalAngle(const eckit::Configuration &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  std::string component_;  // has two valid values: "cos"  and "sin"
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_ORBITALANGLE_H_
