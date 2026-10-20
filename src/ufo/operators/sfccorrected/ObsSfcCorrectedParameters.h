/*
 * (C) Crown Copyright 2024, Met Office
 * (C) Copyright 2024 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
#define UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// enum type for surface correction type, and ParameterTraitsHelper for it
enum class SfcCorrectionType {
  UKMO,   /// UKMO method uses model surface and 2km data.
  WRFDA,  /// WRFDA method uses model surface and level 1 data.
  GSL     /// GSL method uses model surface and level 1 data.
};
struct SfcCorrectionTypeParameterTraitsHelper {
  using EnumType = SfcCorrectionType;
  static constexpr char enumTypeName[] = "SfcCorrectionType";
  static constexpr util::NamedEnumerator<SfcCorrectionType> namedValues[] = {
    { SfcCorrectionType::UKMO, "UKMO" },
    { SfcCorrectionType::WRFDA, "WRFDA" },
    { SfcCorrectionType::GSL, "GSL" }
  };
};

/// Linearization mode for the GSL lapse rate in TL/AD.
/// MatchSimobs: differentiate the configured lapse-rate calculation
/// (consistent linearization).
/// NoAdjustment: treat trajectory-time lapse rate as constant in TL/AD,
/// dropping contribution from model temperature levels used to compute it.
enum class GslTladLapseRateOption {
  MatchSimobs,
  NoAdjustment
};
struct GslTladLapseRateOptionParameterTraitsHelper {
  using EnumType = GslTladLapseRateOption;
  static constexpr char enumTypeName[] = "GslTladLapseRateOption";
  static constexpr util::NamedEnumerator<GslTladLapseRateOption> namedValues[] = {
    { GslTladLapseRateOption::MatchSimobs, "MatchSimobs" },
    { GslTladLapseRateOption::NoAdjustment, "NoAdjustment" }
  };
};


/// Enum type for GSL lapse rate options, only used when SfcCorrectionType is set to GSL
enum class GslLapseRateOption {
  Constant,       /// Use a constant lapse rate value
  Local,          /// Calculate lapse rate locally from model levels
  NoAdjustment    /// No temperature adjustment
};
struct GslLapseRateOptionParameterTraitsHelper {
  using EnumType = GslLapseRateOption;
  static constexpr char enumTypeName[] = "GslLapseRateOption";
  static constexpr util::NamedEnumerator<GslLapseRateOption> namedValues[] = {
    { GslLapseRateOption::Constant, "Constant" },
    { GslLapseRateOption::Local, "Local" },
    { GslLapseRateOption::NoAdjustment, "NoAdjustment" }
  };
};

}  // namespace ufo

namespace oops {

/// Extraction of SfcCorrectionType parameters from config
template <>
struct ParameterTraits<ufo::SfcCorrectionType> :
    public EnumParameterTraits<ufo::SfcCorrectionTypeParameterTraitsHelper>
{};

/// Extraction of GslLapseRateOption parameters from config
template <>
struct ParameterTraits<ufo::GslLapseRateOption> :
    public EnumParameterTraits<ufo::GslLapseRateOptionParameterTraitsHelper>
{};

/// Extraction of GslTladLapseRateOption parameters from config
template <>
struct ParameterTraits<ufo::GslTladLapseRateOption> :
    public EnumParameterTraits<ufo::GslTladLapseRateOptionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/**
 * GSL-specific correction parameters
 * Only used when SfcCorrectionType is set to GSL
 */
class GslCorrectionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GslCorrectionParameters, Parameters)

 public:
  /// Lapse rate calculation method
  oops::Parameter<GslLapseRateOption> temperatureLapseRateOption{
    "temperature lapse rate option",
    "Method to determine lapse rate for surface temperature correction: "
    "('Constant' or 'Local' or 'NoAdjustment')",
    GslLapseRateOption::Local,
    this};

  /// Constant lapse rate value (K/km)
  /// Only used when lapseRateOption = Constant
  oops::Parameter<float> temperatureLapseRateValue{
    "temperature lapse rate",
    "Fixed lapse rate (K/km) used to adjust observed surface temperature "
    "to model surface level. Only used when temperature lapse rate option = Constant."
    "Default: 9.8 K/km (standard adiabatic lapse rate)",
    9.8,
    this
  };

  /// Local lapse rate calculation parameters
  /// Only used when lapseRateOption = Local
  oops::Parameter<int> temperatureLocalLapseRateLevel{
    "temperature local lapse rate level",
    "Highest model level used to calculate local lapse rate."
    "Only used when temperature lapse rate option = Local",
    5,
    this
  };

  /// Apply thresholds to local lapse rate
  /// Only used when lapseRateOption = Local
  oops::Parameter<bool> temperatureLapseRateThreshold{
    "temperature lapse rate threshold",
    "Apply min/max thresholds to calculated local lapse rate."
    "Only used when temperature lapse rate option = Local",
    true,
    this
  };

  /// Minimum threshold for local lapse rate (K/km)
  /// Only used when lapseRateOption = Local and applyThreshold = true
  oops::Parameter<float> minThreshold{
    "min threshold",
    "Minimum lapse rate (K/km) allowed when calculated locally."
    "Only used when temperature lapse rate option = Local and "
    "temperature lapse rate threshold = true",
    0.5,
    this
  };

  /// Maximum threshold for local lapse rate (K/km)
  /// Only used when lapseRateOption = Local and applyThreshold = true
  oops::Parameter<float> maxThreshold{
    "max threshold",
    "Maximum lapse rate (K/km) allowed when calculated locally."
    "Only used when temperature lapse rate option = Local and "
    "temperature lapse rate threshold = true",
    10.0,
    this
  };

  /// Linearization mode for the lapse rate in TL/AD. Defaults to
  /// MatchSimobs (consistent linearization). Set to NoAdjustment to
  /// freeze the lapse rate at trajectory time in TL/AD.
  oops::Parameter<GslTladLapseRateOption> tladLapseRateOption{
    "tlad lapse rate option",
    "Linearization mode for the lapse rate in TL/AD: 'MatchSimobs' "
    "(consistent linearization, default) or 'NoAdjustment' (treat the "
    "trajectory lapse rate as constant in TL/AD.",
    GslTladLapseRateOption::MatchSimobs,
    this
  };

  /// If true, and lapse rate has zero derivative w.r.t. the model state
  /// (temperature lapse rate option = Constant or NoAdjustment, or
  /// tlad lapse rate option = NoAdjustment freezing a Local lapse rate),
  /// add the lowest model level's temperature increment to the
  /// air_temperature_at_2m increment in TL/AD.
  /// Defaults to false, exact linearization (dhofx = dT_2m).
  /// If the model's TL/AD does not carry an independent
  /// active increment for air_temperature_at_2m, this can be used
  /// so the observation can reach the 3D prognostic field.
  oops::Parameter<bool> tladSpreadToLowestLevel{
    "tlad spread to lowest level",
    "When the lapse rate has zero derivative w.r.t. the model state "
    "(Constant, NoAdjustment, or tlad lapse rate option = NoAdjustment), "
    "add the lowest model level's air_temperature increment to the "
    "air_temperature_at_2m increment in TL/AD, so the observation "
    "increment still reaches the 3D prognostic field",
    false,
    this
  };
};

/// Configuration options recognized by the SfcCorrected operator.
class ObsSfcCorrectedParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSfcCorrectedParameters, ObsOperatorParametersBase)

 public:
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
      "variables",
      "List of variables to be simulated which must be a subset of the simulated variables "
      "in the ObsSace",
      this};

  oops::Parameter<std::string> geovarGeomZ{
      "geovar_geomz",
      "Model variable for height of vertical levels, geopotential heights will be converted",
      "height_above_mean_sea_level",
      this};

  oops::Parameter<std::string> geovarSfcGeomZ{
      "geovar_sfc_geomz",
      "Model variable for surface height, geopotential heights will be converted",
      "height_above_mean_sea_level_at_surface",
      this};

  oops::Parameter<std::string> obsHeightName{
      "station_altitude",
      "stationElevation",
      this};

  oops::Parameter<ufo::Variable> obsTemperature{
      "observed temperature variable",
      "Observed air temperature variable (group/variable format)",
      ufo::Variable("ObsValue/airTemperatureAt2M"),
      this};

  oops::Parameter<ufo::Variable> obsSpecificHumidity{
      "observed specific humidity variable",
      "Observed specific humidity variable (group/variable format)",
      ufo::Variable("ObsValue/specificHumidityAt2M"),
      this};

  oops::Parameter<SfcCorrectionType> correctionType{
      "correction scheme to use",
      "Scheme used for correction ('WRFDA' or 'UKMO' or 'GSL')",
      SfcCorrectionType::WRFDA,
      this};

  /// GSL-specific configuration parameters
  /// Only used when correctionType = GSL
  oops::OptionalParameter<GslCorrectionParameters> gslParams{
    "gsl parameters",
    "GSL-specific surface correction parameters."
    "Only used when correction_scheme = GSL",
    this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
