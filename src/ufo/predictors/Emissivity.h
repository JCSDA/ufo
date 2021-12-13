/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_EMISSIVITY_H_
#define UFO_PREDICTORS_EMISSIVITY_H_

#include "ufo/predictors/PredictorBase.h"

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class Emissivity : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef EmptyPredictorParameters Parameters_;

  Emissivity(const Parameters_ &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_EMISSIVITY_H_
