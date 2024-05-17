/*
 * (C) Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_READBIAS_H_
#define UFO_PREDICTORS_READBIAS_H_

#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
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

/// Configuration parameters of the ReadBias predictor.
class ReadBiasParameters: public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ReadBiasParameters, PredictorParametersBase);

 public:
  /// Name of the group containing the values to be read.
  oops::Parameter<std::string> groupName{"group name", "ObsBias", this};
};

class ReadBias : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef ReadBiasParameters Parameters_;

  ReadBias(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;
 private:
  std::string group_name_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_READBIAS_H_
