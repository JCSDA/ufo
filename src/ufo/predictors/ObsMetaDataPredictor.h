/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_OBSMETADATAPREDICTOR_H_
#define UFO_PREDICTORS_OBSMETADATAPREDICTOR_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace oops {
  class ObsVariables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the predictor.
class ObsMetaDataPredictorParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsMetaDataPredictorParameters, PredictorParametersBase);

 public:
  /// Power to which to raise the predictor. By default, 1.
  ///
  /// \note If this option is set, a suffix containing its value (even if it's equal to 1) will be
  /// appended to the predictor name.
  oops::OptionalParameter<int> order{"order", this};

  /// Name of the variable (from the obs MetaData group) containing the predictor.
  oops::RequiredParameter<std::string> varName{"variable", this};
};

// -----------------------------------------------------------------------------

class ObsMetaDataPredictor : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef ObsMetaDataPredictorParameters Parameters_;

  ObsMetaDataPredictor(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  std::string variable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_OBSMETADATAPREDICTOR_H_
