/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_

#include <string>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Parameters for the observed variable configuration
class ObservedVariableParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObservedVariableParameters, Parameters)

 public:
  /// Full variable name (group/variable format)
  oops::Parameter<std::string> name{"name", "ObsValue/airTemperatureAt2M", this};
};

/// \brief Parameters for the model height configuration
class ModelHeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightParameters, Parameters)

 public:
  /// From: "terrain" (model orography) or "lowestModelLevel"
  oops::Parameter<std::string> from{"from", "terrain", this};

  /// Offset added to the model height in meters
  oops::Parameter<float> offsetInMeters{"offset in meters", 0.0f, this};
};

/// \brief Parameters for lapse rate bounds
class LapseRateBoundsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LapseRateBoundsParameters, Parameters)

 public:
  /// Minimum allowed lapse rate in K/km
  oops::RequiredParameter<float> minvalue{"minvalue", this};
  /// Maximum allowed lapse rate in K/km
  oops::RequiredParameter<float> maxvalue{"maxvalue", this};
};

/// \brief Parameters for the temperature lapse rate configuration
class TemperatureLapseRateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TemperatureLapseRateParameters, Parameters)

 public:
  /// Lapse rate scheme: "constant" or "local"
  oops::Parameter<std::string> scheme{"scheme", "constant", this};

  /// Constant lapse rate value in K/km (used when option is "constant")
  oops::Parameter<float> value{"value", 6.5f, this};

  /// Model level index1 (counted from surface, 1-based)
  oops::Parameter<int> level1{"level1", 1, this};

  /// Model level index2 (counted from surface, 1-based)
  oops::Parameter<int> level2{"level2", 5, this};

  /// Optional bounds to clamp the local lapse rate (K/km)
  oops::OptionalParameter<LapseRateBoundsParameters> bounds{"bounds", this};
};

/// \brief Options controlling ModelHeightAdjustedAirTemperature ObsFunction
class ModelHeightAdjustedAirTemperatureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightAdjustedAirTemperatureParameters, Parameters)

 public:
  /// Input observation station height to be used
  oops::RequiredParameter<Variable> elevation{"elevation", this};

  /// Observed variable configuration (optional; defaults to ObsValue/airTemperatureAt2M)
  oops::OptionalParameter<ObservedVariableParameters> observedVariable{
      "observed variable", this};

  /// Model height configuration (optional; defaults to terrain with 0m offset)
  oops::OptionalParameter<ModelHeightParameters> modelHeight{
      "model height", this};

  /// Temperature lapse rate configuration (optional; defaults to constant 6.5 K/km)
  oops::OptionalParameter<TemperatureLapseRateParameters> lapseRate{
      "temperature lapse rate in K/km", this};
};

/// \brief  Function to adjust a surface temperature observation value from station height
/// to model surface height using a lapse rate correction.
///
/// Supports two modes:
/// - "constant" (default): Uses a configurable lapse rate (default 6.5 K/km)
/// - "local": Computes lapse rate from model temperature profile between two levels
///
/// The observation variable defaults to airTemperatureAt2M but is configurable.

class ModelHeightAdjustedAirTemperature : public ObsFunctionBase<float> {
 public:
    explicit ModelHeightAdjustedAirTemperature(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ModelHeightAdjustedAirTemperatureParameters parameters_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_
