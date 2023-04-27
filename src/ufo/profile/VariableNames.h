/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_VARIABLENAMES_H_
#define UFO_PROFILE_VARIABLENAMES_H_

namespace ufo
{

struct VariableNames
{
  // Observation values

  static constexpr const char* const obs_air_pressure = "MetaData/pressure";
  static constexpr const char* const obs_air_temperature = "ObsValue/airTemperature";
  static constexpr const char* const obs_relative_humidity = "ObsValue/relativeHumidity";
  static constexpr const char* const obs_eastward_wind = "ObsValue/windEastward";
  static constexpr const char* const obs_northward_wind = "ObsValue/windNorthward";
  static constexpr const char* const obs_geopotential_height = "ObsValue/height";
  static constexpr const char* const obs_dew_point_temperature = "ObsValue/dewPointTemperature";

  // Derived observation values

  static constexpr const char* const obs_derived_air_pressure = "DerivedMetaData/pressure";

  // Observation errors

  static constexpr const char* const obserr_air_temperature = "ObsErrorData/airTemperature";
  static constexpr const char* const obserr_relative_humidity = "ObsErrorData/relativeHumidity";
  static constexpr const char* const obserr_eastward_wind = "ObsErrorData/windEastward";
  static constexpr const char* const obserr_northward_wind = "ObsErrorData/windNorthward";
  static constexpr const char* const obserr_geopotential_height = "ObsErrorData/height";
  static constexpr const char* const obserr_dew_point_temperature =
    "ObsErrorData/dewPointTemperature";

  // HofX

  static constexpr const char* const hofx_air_temperature = "HofX/airTemperature";
  static constexpr const char* const hofx_relative_humidity = "HofX/relativeHumidity";
  static constexpr const char* const hofx_eastward_wind = "HofX/windEastward";
  static constexpr const char* const hofx_northward_wind = "HofX/windNorthward";
  static constexpr const char* const hofx_geopotential_height = "HofX/height";
  static constexpr const char* const hofx_dew_point_temperature = "HofX/dewPointTemperature";

  // Probability of gross error

  static constexpr const char* const pge_air_temperature = "GrossErrorProbability/airTemperature";
  static constexpr const char* const pge_relative_humidity =
    "GrossErrorProbability/relativeHumidity";
  static constexpr const char* const pge_eastward_wind = "GrossErrorProbability/windEastward";
  static constexpr const char* const pge_northward_wind = "GrossErrorProbability/windNorthward";
  static constexpr const char* const pge_geopotential_height =
    "GrossErrorProbability/height";

  // MetaData

  static constexpr const char* const station_ID = "MetaData/stationIdentification";
  static constexpr const char* const ObsType = "MetaData/observationSubTypeNum";
  static constexpr const char* const Latitude = "MetaData/latitude";
  static constexpr const char* const Longitude = "MetaData/longitude";
  static constexpr const char* const Zstation = "MetaData/stationElevation";
  static constexpr const char* const LevelType = "MetaData/levelType";
  static constexpr const char* const InstrType = "MetaData/instrumentIdentifier";
  static constexpr const char* const extended_obs_space = "MetaData/extendedObsSpace";

  // QC flags

  static constexpr const char* const qcflags_observation_report = "QCFlags/observationReport";
  static constexpr const char* const qcflags_air_temperature = "QCFlags/airTemperature";
  static constexpr const char* const qcflags_relative_humidity = "QCFlags/relativeHumidity";
  static constexpr const char* const qcflags_geopotential_height = "QCFlags/height";
  static constexpr const char* const qcflags_eastward_wind = "QCFlags/windEastward";
  static constexpr const char* const qcflags_northward_wind = "QCFlags/windNorthward";
  static constexpr const char* const qcflags_wind_profiler = "QCFlags/windProfiler";

  // Diagnostic flags

  static constexpr const char* const diagflags_profile_interpolation_eastward_wind =
    "DiagnosticFlags/Profile/Interpolation/windEastward";
  static constexpr const char* const diagflags_profile_standard_level_eastward_wind =
    "DiagnosticFlags/Profile/StandardLevel/windEastward";

  // Counters

  static constexpr const char* const counter_NumAnyErrors = "Counters/NumAnyErrors";
  static constexpr const char* const counter_NumSamePErrObs = "Counters/NumSamePErrObs";
  static constexpr const char* const counter_NumSuperadiabat = "Counters/NumSuperadiabat";
  static constexpr const char* const counter_Num925Miss = "Counters/Num925Miss";
  static constexpr const char* const counter_Num100Miss = "Counters/Num100Miss";
  static constexpr const char* const counter_NumStdMiss = "Counters/NumStdMiss";
  static constexpr const char* const counter_NumHydErrObs = "Counters/NumHydErrObs";
  static constexpr const char* const counter_NumIntHydErrors = "Counters/NumIntHydErrors";
  static constexpr const char* const counter_NumInterpErrors = "Counters/NumInterpErrors";
  static constexpr const char* const counter_NumInterpErrObs = "Counters/NumInterpErrObs";
  static constexpr const char* const counter_NumSignChange = "Counters/NumSignChange";
  static constexpr const char* const counter_TotCProfs = "Counters/TotCProfs";
  static constexpr const char* const counter_TotHProfs = "Counters/TotHProfs";
  static constexpr const char* const counter_TotCFlags = "Counters/TotCFlags";
  static constexpr const char* const counter_TotHFlags = "Counters/TotHFlags";
  static constexpr const char* const counter_TotLFlags = "Counters/TotLFlags";
  static constexpr const char* const counter_NumGapsT = "Counters/NumGapsT";
  static constexpr const char* const counter_NumGapsU = "Counters/NumGapsU";
  static constexpr const char* const counter_NumGapsUWP = "Counters/NumGapsUWP";
  static constexpr const char* const counter_NumGapsRH = "Counters/NumGapsRH";

  // Corrections

  static constexpr const char* const obscorrection_air_temperature = "Corrections/airTemperature";
  static constexpr const char* const obscorrection_geopotential_height =
    "Corrections/height";

  // Intermediate values

  static constexpr const char* const DC = "MetaData/DC";
  static constexpr const char* const ETol = "MetaData/ETol";
  static constexpr const char* const D = "MetaData/D";
  static constexpr const char* const E = "MetaData/E";
  static constexpr const char* const HydError = "MetaData/HydError";
  static constexpr const char* const PBottom = "MetaData/PBottom";
  static constexpr const char* const StdLev = "MetaData/StdLev";
  static constexpr const char* const SigAbove = "MetaData/SigAbove";
  static constexpr const char* const SigBelow = "MetaData/SigBelow";
  static constexpr const char* const IndStd = "MetaData/IndStd";
  static constexpr const char* const LevErrors = "MetaData/LevErrors";
  static constexpr const char* const tInterp = "MetaData/tInterp";
  static constexpr const char* const uInterp = "MetaData/uInterp";
  static constexpr const char* const vInterp = "MetaData/vInterp";
  static constexpr const char* const LogP = "MetaData/LogP";
  static constexpr const char* const NumStd = "MetaData/NumStd";
  static constexpr const char* const NumSig = "MetaData/NumSig";
  static constexpr const char* const Press = "MetaData/Press";
  static constexpr const char* const Temp = "MetaData/Temp";
  static constexpr const char* const rh = "MetaData/rh";
  static constexpr const char* const td = "MetaData/td";
  static constexpr const char* const tbk = "MetaData/tbk";
  static constexpr const char* const rhbk = "MetaData/rhbk";
  static constexpr const char* const FlagH = "MetaData/FlagH";
  static constexpr const char* const Indx = "MetaData/Indx";

  // GeoVaLs

  static constexpr const char* const geovals_orog = "surface_altitude";
  static constexpr const char* const geovals_pressure = "air_pressure";
  static constexpr const char* const geovals_pressure_rho = "air_pressure_levels";
  static constexpr const char* const geovals_pressure_rho_minus_one =
    "air_pressure_levels_minus_one";
  static constexpr const char* const geovals_height = "height";
  static constexpr const char* const geovals_height_rho = "height_levels";
  static constexpr const char* const geovals_height_rho_minus_one = "height_levels_minus_one";
  static constexpr const char* const geovals_potential_temperature = "potential_temperature";
  static constexpr const char* const geovals_air_temperature = "air_temperature";
  static constexpr const char* const geovals_surface_pressure = "surface_pressure";
  static constexpr const char* const geovals_relative_humidity = "relative_humidity";

  // GeoVaLs used in validation

  static constexpr const char* const geovals_testreference_logP = "LogPB";
  static constexpr const char* const geovals_testreference_ExnerP = "ExnerPB";
  static constexpr const char* const geovals_testreference_logP_rho = "LogPA";
  static constexpr const char* const geovals_testreference_ExnerP_rho = "ExnerPA";
  static constexpr const char* const geovals_testreference_air_temperature =
    "average_air_temperature";
  static constexpr const char* const geovals_testreference_eastward_wind =
    "average_eastward_wind";
  static constexpr const char* const geovals_testreference_northward_wind =
    "average_northward_wind";
  static constexpr const char* const geovals_testreference_relative_humidity =
    "average_relative_humidity";
  static constexpr const char* const geovals_testreference_air_temperature_qcflags =
    "average_air_temperature_flags";
  static constexpr const char* const geovals_testreference_eastward_wind_qcflags =
    "average_eastward_wind_flags";
  static constexpr const char* const geovals_testreference_northward_wind_qcflags =
    "average_northward_wind_flags";
  static constexpr const char* const geovals_testreference_relative_humidity_qcflags =
    "average_relative_humidity_flags";

  // Averaged values on model levels

  static constexpr const char* const air_temperature_derived =
    "DerivedObsValue/airTemperature";
  static constexpr const char* const eastward_wind_derived =
    "DerivedObsValue/windEastward";
  static constexpr const char* const northward_wind_derived =
    "DerivedObsValue/windNorthward";
  static constexpr const char* const relative_humidity_derived =
    "DerivedObsValue/relativeHumidity";

  // Derived observation values (used in averaging)

  static constexpr const char* const LogP_derived = "DerivedMetaData/logP";
  static constexpr const char* const bigPgaps_derived = "DerivedMetaData/bigPgaps";

  // Derived model values (used in averaging)

  static constexpr const char* const modellevels_logP_derived = "DerivedModelValue/LogPB";
  static constexpr const char* const modellevels_ExnerP_derived = "DerivedModelValue/ExnerPB";
  static constexpr const char* const modellevels_air_temperature_derived =
    "DerivedModelValue/airTemperature";

  // Derived model values on rho levels (used in averaging)

  static constexpr const char* const modellevels_logP_rho_derived =
    "DerivedModelValue/LogPA";
  static constexpr const char* const modellevels_logPWB_rho_derived =
    "DerivedModelValue/LogPWB";
  static constexpr const char* const modellevels_ExnerP_rho_derived =
    "DerivedModelValue/ExnerPA";
};

}  // namespace ufo

#endif  // UFO_PROFILE_VARIABLENAMES_H_
