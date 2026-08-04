/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_STABLELAYERSCLOUDTOPPRESSURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_STABLELAYERSCLOUDTOPPRESSURE_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling StableLayersCloudTopPressure ObsFunction
class StableLayersCloudTopPressureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StableLayersCloudTopPressureParameters, Parameters)

 public:
  /// Set of active channels.
  oops::RequiredParameter<std::string> channels{"channels", this};

  /// Conditions used to select locations at which the StableLayersCloudTopPressure
  /// ObsFunction should be applied.
  /// If not specified, all locations will be selected.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// Tropopause level.
  oops::RequiredParameter<std::string> tropopauseLevel
    {"tropopause level",
     "Value of tropopause level.",
     this};

  /// Bias corrected brightness temperature for window channel.
  /// Usually of the form OneDVar/brightnessTemperature_9 for window channel 9.
  oops::RequiredParameter<std::string> correctedBrightnessTemperature
    {"corrected brightness temperature",
     "Name of corrected brightness temperature.",
     this};

  /// Window channel used for selecting RTTOV brightness temperature.
  /// This channel should be the same as the channel passed into correctedBrightnessTemperature.
  oops::Parameter<int> windowChannel
    {"window channel",
     "Integer value of window channel (default=9).",
     9,
     this};

  /// Apply brightness temperature constraint.
  oops::Parameter<bool> useBTConstraint
    {"brightness temperature constraint",
     "Constraint that determines whether brightness temperature is used to detect "
     "stable layers and to set the stable weights (default=false).",
     false,
     this};

  /// Warm temperature limit.
  oops::Parameter<float> tempLimitWarm
    {"temperature limit warm",
     "Limit for the positive difference between the bias corrected brightness temperature and "
     "the brightness temperature of the level (default=1.0).",
     1.0,
     this};

  /// Cold temperature limit.
  oops::Parameter<float> tempLimitCold
    {"temperature limit cold",
     "Limit for the negative difference between the bias corrected brightness temperature and "
     "the brightness temperature of the level (default=-1.0).",
     -1.0,
     this};

  /// Stable density.
  oops::Parameter<float> stableDensity
    {"stable density",
     "Parameter that controls the density of the lapse rate weights (default=1.0).",
     1.0,
     this};

  /// Relative humidity density, in units of model relative humidity (i.e. as a
  /// fraction rather than a percentage).
  oops::Parameter<float> relativeHumidityDensity{
      "relative humidity density as a fraction",
      "Parameter that controls the density of the relative humidity weights, "
      "in units of model relative humidity (i.e. as a fraction rather than a "
      "percentage) (default=1.0).",
      1.0, this};

  /// Relative humidity offset.
  oops::Parameter<float> relativeHumidityOffset
    {"relative humidity offset",
     "Parameter that controls the offset of the relative humidity weights (default=0.0).",
     0.0,
     this};

  /// Relative humidity minimum.
  oops::Parameter<float> relativeHumidityMinimum
    {"relative humidity minimum",
     "Parameter that controls the minimum of the relative humidity weights (default=0.0).",
     0.0,
     this,
     {oops::minConstraint(0.0f)}};
};

// -----------------------------------------------------------------------------

/// \brief Calculates the cloud top pressure, brightness temperature at the cloud top and the
/// weighted standard deviation about the cloud top pressure using the Stable Layers method.
///
/// There are five steps to the Stable Layers method. First, each layer of the atmosphere undergoes
/// up to three tests to see if it is a possible stable layer. The first test checks if the lapse
/// rate is less than the saturated adiabatic lapse rate. The second and third tests, only
/// conducted if useBTConstraint is set to true, check if the bias corrected brightness temperature
/// is within a certain range of the brightness temperature at the level. The parameters
/// tempLimitWarm and tempLimitCold are used to set a tolerance for the second and third tests.
/// Then, each layer that passes the tests is assigned a stability weight. The parameters
/// stableDensity, relativeHumidityDensity, relativeHumidityOffset and relativeHumidityMinimum
/// are used to control the calculation of the stability weight.
/// Next, the cloud top pressure is calculated using a quadratic fit to the stability weights.
/// The brightness temperature at the cloud top is then found by interpolating the model
/// brightness temperature at the cloud top pressure.
/// Finally, the weighted standard deviation about the cloud top pressure is calculated.
class StableLayersCloudTopPressure : public ObsFunctionBase<float> {
 public:
  explicit StableLayersCloudTopPressure(const eckit::LocalConfiguration &);
  ~StableLayersCloudTopPressure();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  StableLayersCloudTopPressureParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;

  double computeSALR(const double, const double) const;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_STABLELAYERSCLOUDTOPPRESSURE_H_
