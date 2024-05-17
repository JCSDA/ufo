/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_THICKNESS_H_
#define UFO_PREDICTORS_THICKNESS_H_

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the thickness predictor.
class ThicknessParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ThicknessParameters, PredictorParametersBase);

 public:
  /// Pressure value (Pa) at the top of the required thickness layer
  oops::RequiredParameter<double> layerTop{"layer top", this};

  /// Pressure value (Pa) at the bottom of the required thickness layer
  oops::RequiredParameter<double> layerBase{"layer base", this};

  /// Climatological mean of predictor
  oops::RequiredParameter<double> mean{"mean", this};

  /// Climatological standard deviation of predictor
  oops::RequiredParameter<double> stDev{"standard deviation", this};
};

// -----------------------------------------------------------------------------

/**
 *This thickness predictor calculates the thickness of a specified pressure level interval.
 *The thickness is calculated as the difference between the geopotential heights at two pressure
 *levels. This requires the integration of the temperature with respect to log pressure.
 *The thicknesses are based on geovals temperature and pressure profiles
 */


class Thickness : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef ThicknessParameters Parameters_;

  Thickness(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_THICKNESS_H_
