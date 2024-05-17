/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASPARAMETERS_H_
#define UFO_OBSBIASPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsObsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace ufo {

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// PredictorParametersBase.
class PredictorParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PredictorParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of PredictorParametersBase configuring
  /// a predictor. The type of the subclass is determined by the value of the "name" key in the
  /// Configuration object from which this object is deserialized.
  oops::RequiredPolymorphicParameter<PredictorParametersBase, PredictorFactory>
    predictorParameters{"name", this};
};

class StaticOrVariationalBCParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StaticOrVariationalBCParameters, Parameters)

 public:
  /// Each element of this list is used to configure a separate predictor.
  oops::Parameter<std::vector<PredictorParametersWrapper>> predictors{"predictors", {}, this};
};

class ObsBiasCovariancePriorInflationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasCovariancePriorInflationParameters, Parameters)

 public:
  oops::RequiredParameter<double> ratio{"ratio", this};
  oops::RequiredParameter<double> ratioForSmallDataset{"ratio for small dataset", this};
};

class ObsBiasCovariancePriorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasCovariancePriorParameters, Parameters)

 public:
  oops::RequiredParameter<ObsBiasCovariancePriorInflationParameters> inflation{"inflation", this};
  oops::OptionalParameter<std::string> inputFile{"input file", this};
};

class ObsBiasCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasCovarianceParameters, Parameters)

 public:
  /// Default smallest variance value
  static double defaultSmallestVariance() { return 1.0e-6; }
  /// Default largest variance value
  static double defaultLargestVariance() { return 10.0; }
  /// Default largest analysis error variance
  static double defaultLargestAnalysisVariance() { return 10000.0; }
  /// Default step size
  static double defaultStepSize() { return 1.e-4; }

  oops::RequiredParameter<size_t> minimalRequiredObsNumber{
    "minimal required obs number", this};
  oops::Parameter<std::vector<double>> varianceRange{
    "variance range", {defaultSmallestVariance(), defaultLargestVariance()}, this};
  oops::Parameter<double> stepSize{
    "step size", defaultStepSize(), this};
  oops::Parameter<double> largestAnalysisVariance{
    "largest analysis variance", defaultLargestAnalysisVariance(), this};

  oops::OptionalParameter<ObsBiasCovariancePriorParameters> prior{
    "prior", this};

  oops::OptionalParameter<std::string> outputFile{"output file", this};
};

/// Parameters influencing the bias correction process.
class ObsBiasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasParameters, Parameters)

 public:
  /// List of predictors with unit coefficients (unaffected by VarBC).
  oops::Parameter<StaticOrVariationalBCParameters> staticBC{"static bc", {}, this};
  /// List of predictors with coefficients determined by variational analysis (VarBC).
  oops::Parameter<StaticOrVariationalBCParameters> variationalBC{"variational bc", {}, this};
  /// Boolean determining whether records are bias-corrected independently from each other.
  oops::Parameter<bool> BiasCorrectionByRecord{"bc by record", false, this};
  /// List of variables (and channels) that is not bias corrected
  oops::Parameter<oops::ObsVariables> variablesNoBC{"variables without bc", {}, this};
  /// Path to a NetCDF file containing initial values of the coefficients of predictors used
  /// in VarBC.
  oops::OptionalParameter<std::string> inputFile{"input file", this};
  /// Path to a NetCDF file containing final values of the coefficients of predictors used
  /// in VarBC.
  oops::OptionalParameter<std::string> outputFile{"output file", this};
  /// Path to a NetCDF file containing increments of the coefficients of predictors used
  /// in VarBC.
  oops::OptionalParameter<std::string> outputFileInc{"increment output file", this};
  /// Options controlling the covariance matrix.
  oops::OptionalParameter<ObsBiasCovarianceParameters> covariance{"covariance", this};
};

}  // namespace ufo

#endif  // UFO_OBSBIASPARAMETERS_H_
