/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_LAPSERATE_H_
#define UFO_PREDICTORS_LAPSERATE_H_

#include <map>
#include <string>

#include "oops/util/parameters/OptionalParameter.h"
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

/// Configuration parameters of the LapseRate predictor.
class LapseRateParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(LapseRateParameters, PredictorParametersBase);

 public:
  /// Path to an input file.
  oops::RequiredParameter<std::string> tlapse{"tlapse", this};

  /// Power to which to raise the lapse rate. By default, 1.
  ///
  /// \note If this option is set, a suffix containing its value (even if it's equal to 1) will be
  /// appended to the predictor name.
  oops::OptionalParameter<int> order{"order", this};
};

// -----------------------------------------------------------------------------

class LapseRate : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef LapseRateParameters Parameters_;

  LapseRate(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  std::map<int, float> tlapmean_;  // <channel, tlaps>
  int order_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_LAPSERATE_H_
