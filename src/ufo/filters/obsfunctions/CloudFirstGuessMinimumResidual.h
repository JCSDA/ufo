/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLOUDFIRSTGUESSMINIMUMRESIDUAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLOUDFIRSTGUESSMINIMUMRESIDUAL_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Parameters for the Cloud Top Minimum Residual
///
class CloudFirstGuessMinimumResidualParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CloudFirstGuessMinimumResidualParameters, Parameters)

 public:
  /// Set of channels used in the calculation of the cost function
  oops::RequiredParameter<std::string> channels{"channels", this};

  /// Output directory
  oops::Parameter<std::string> outputGroup{"output group", "MetaData", this};

  /// Output name for cloud top pressure
  oops::Parameter<std::string> cloudTopPressureName{"output name for cloud top pressure",
                                                    "pressureAtTopOfCloud", this};

  /// Output name for cloud fraction
  oops::Parameter<std::string> cloudFractionName{"output name for cloud fraction",
                                                    "cloudAmount", this};

  /// Minimum pressure for the cloud in Pa
  oops::Parameter<float> minCloudPressure{"minimum cloud top pressure", 10000.0f, this};

  /// The group in which the ObsBias is located - this allows for it to be provided from
  /// the ObsSpace for testing
  oops::Parameter<std::string> obsBiasGroup{"obs bias group", "ObsBiasData", this};
};

///
/// \brief Cloud Top using Minimum Residual method
///
class CloudFirstGuessMinimumResidual : public ObsFunctionBase<float> {
 public:
  explicit CloudFirstGuessMinimumResidual(const eckit::LocalConfiguration &);
  ~CloudFirstGuessMinimumResidual();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return invars_;}
 private:
  CloudFirstGuessMinimumResidualParameters options_;
  std::vector<int> channels_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLOUDFIRSTGUESSMINIMUMRESIDUAL_H_
