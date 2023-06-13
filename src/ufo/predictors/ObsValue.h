/*
 * (C) Copyright 2023 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_OBSVALUE_H_
#define UFO_PREDICTORS_OBSVALUE_H_

#include <string>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/predictors/PredictorBase.h"

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the ObsValue predictor.
class ObsValueParameters: public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsValueParameters, PredictorParametersBase);

 public:
  /// Name of the variable (from the ObsValue group) containing the observation value.
  oops::RequiredParameter<std::string> varName{"variable name", this};
};

class ObsValue : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef ObsValueParameters Parameters_;

  ObsValue(const Parameters_ &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;
 private:
  std::string var_name_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_OBSVALUE_H_
