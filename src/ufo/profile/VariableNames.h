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

  static constexpr const char* const obs_air_pressure = "air_pressure@MetaData";
  static constexpr const char* const obs_air_temperature = "air_temperature@ObsValue";
  static constexpr const char* const obs_relative_humidity = "relative_humidity@ObsValue";
  static constexpr const char* const obs_eastward_wind = "eastward_wind@ObsValue";
  static constexpr const char* const obs_northward_wind = "northward_wind@ObsValue";
  static constexpr const char* const obs_geopotential_height = "height@ObsValue";
  static constexpr const char* const obs_dew_point_temperature = "dew_point_temperature@ObsValue";

  // Observation errors

  static constexpr const char* const obserr_air_temperature = "air_temperature@ObsErrorData";
  static constexpr const char* const obserr_relative_humidity = "relative_humidity@ObsErrorData";
  static constexpr const char* const obserr_eastward_wind = "eastward_wind@ObsErrorData";
  static constexpr const char* const obserr_northward_wind = "northward_wind@ObsErrorData";
  static constexpr const char* const obserr_geopotential_height = "height@ObsErrorData";
  static constexpr const char* const obserr_dew_point_temperature =
    "dew_point_temperature@ObsErrorData";

  // HofX

  static constexpr const char* const hofx_air_temperature = "air_temperature@HofX";
  static constexpr const char* const hofx_relative_humidity = "relative_humidity@HofX";
  static constexpr const char* const hofx_eastward_wind = "eastward_wind@HofX";
  static constexpr const char* const hofx_northward_wind = "northward_wind@HofX";
  static constexpr const char* const hofx_geopotential_height = "height@HofX";
  static constexpr const char* const hofx_dew_point_temperature = "dew_point_temperature@HofX";

  // Background errors

  static constexpr const char* const bkgerr_air_temperature =
    "air_temperature_background_error@ObsDiag";
  static constexpr const char* const bkgerr_relative_humidity =
    "relative_humidity_background_error@ObsDiag";
  static constexpr const char* const bkgerr_eastward_wind =
    "eastward_wind_background_error@ObsDiag";
  static constexpr const char* const bkgerr_northward_wind =
    "northward_wind_background_error@ObsDiag";
  static constexpr const char* const bkgerr_geopotential_height =
    "height_background_error@ObsDiag";
  static constexpr const char* const bkgerr_dew_point_temperature =
    "dew_point_temperature_background_error@ObsDiag";

  // Probability of gross error

  static constexpr const char* const pge_air_temperature = "air_temperature@GrossErrorProbability";
  static constexpr const char* const pge_relative_humidity =
    "relative_humidity@GrossErrorProbability";
  static constexpr const char* const pge_eastward_wind = "eastward_wind@GrossErrorProbability";
  static constexpr const char* const pge_northward_wind = "northward_wind@GrossErrorProbability";
  static constexpr const char* const pge_geopotential_height =
    "height@GrossErrorProbability";

  // MetaData

  static constexpr const char* const station_ID = "station_id@MetaData";
  static constexpr const char* const ObsType = "ops_subtype@MetaData";
  static constexpr const char* const Latitude = "latitude@MetaData";
  static constexpr const char* const Longitude = "longitude@MetaData";
  static constexpr const char* const Zstation = "station_altitude@MetaData";
  static constexpr const char* const LevelType = "level_type@MetaData";
  static constexpr const char* const InstrType = "instrument_type@MetaData";
  static constexpr const char* const extended_obs_space = "extended_obs_space@MetaData";

  // QC flags

  static constexpr const char* const qcflags_observation_report = "observation_report@QCFlags";
  static constexpr const char* const qcflags_air_temperature = "air_temperature@QCFlags";
  static constexpr const char* const qcflags_relative_humidity = "relative_humidity@QCFlags";
  static constexpr const char* const qcflags_geopotential_height = "height@QCFlags";
  static constexpr const char* const qcflags_eastward_wind = "eastward_wind@QCFlags";
  static constexpr const char* const qcflags_northward_wind = "northward_wind@QCFlags";
  static constexpr const char* const qcflags_wind_profiler = "wind_profiler@QCFlags";

  // Diagnostic flags

  static constexpr const char* const diagflags_profile_interpolation_eastward_wind =
    "DiagnosticFlags/Profile/Interpolation/eastward_wind";
  static constexpr const char* const diagflags_profile_standard_level_eastward_wind =
    "DiagnosticFlags/Profile/StandardLevel/eastward_wind";

  // Counters

  static constexpr const char* const counter_NumAnyErrors = "NumAnyErrors@Counters";
  static constexpr const char* const counter_NumSamePErrObs = "NumSamePErrObs@Counters";
  static constexpr const char* const counter_NumSuperadiabat = "NumSuperadiabat@Counters";
  static constexpr const char* const counter_Num925Miss = "Num925Miss@Counters";
  static constexpr const char* const counter_Num100Miss = "Num100Miss@Counters";
  static constexpr const char* const counter_NumStdMiss = "NumStdMiss@Counters";
  static constexpr const char* const counter_NumHydErrObs = "NumHydErrObs@Counters";
  static constexpr const char* const counter_NumIntHydErrors = "NumIntHydErrors@Counters";
  static constexpr const char* const counter_NumInterpErrors = "NumInterpErrors@Counters";
  static constexpr const char* const counter_NumInterpErrObs = "NumInterpErrObs@Counters";
  static constexpr const char* const counter_NumSignChange = "NumSignChange@Counters";
  static constexpr const char* const counter_TotCProfs = "TotCProfs@Counters";
  static constexpr const char* const counter_TotHProfs = "TotHProfs@Counters";
  static constexpr const char* const counter_TotCFlags = "TotCFlags@Counters";
  static constexpr const char* const counter_TotHFlags = "TotHFlags@Counters";
  static constexpr const char* const counter_TotLFlags = "TotLFlags@Counters";
  static constexpr const char* const counter_NumGapsT = "NumGapsT@Counters";
  static constexpr const char* const counter_NumGapsU = "NumGapsU@Counters";
  static constexpr const char* const counter_NumGapsUWP = "NumGapsUWP@Counters";
  static constexpr const char* const counter_NumGapsRH = "NumGapsRH@Counters";

  // Corrections

  static constexpr const char* const obscorrection_air_temperature = "air_temperature@Corrections";
  static constexpr const char* const obscorrection_geopotential_height =
    "height@Corrections";

  // Intermediate values

  static constexpr const char* const DC = "DC@MetaData";
  static constexpr const char* const ETol = "ETol@MetaData";
  static constexpr const char* const D = "D@MetaData";
  static constexpr const char* const E = "E@MetaData";
  static constexpr const char* const HydError = "HydError@MetaData";
  static constexpr const char* const PBottom = "PBottom@MetaData";
  static constexpr const char* const StdLev = "StdLev@MetaData";
  static constexpr const char* const SigAbove = "SigAbove@MetaData";
  static constexpr const char* const SigBelow = "SigBelow@MetaData";
  static constexpr const char* const IndStd = "IndStd@MetaData";
  static constexpr const char* const LevErrors = "LevErrors@MetaData";
  static constexpr const char* const tInterp = "tInterp@MetaData";
  static constexpr const char* const uInterp = "uInterp@MetaData";
  static constexpr const char* const vInterp = "vInterp@MetaData";
  static constexpr const char* const LogP = "LogP@MetaData";
  static constexpr const char* const NumStd = "NumStd@MetaData";
  static constexpr const char* const NumSig = "NumSig@MetaData";
  static constexpr const char* const Press = "Press@MetaData";
  static constexpr const char* const Temp = "Temp@MetaData";
  static constexpr const char* const rh = "rh@MetaData";
  static constexpr const char* const td = "td@MetaData";
  static constexpr const char* const tbk = "tbk@MetaData";
  static constexpr const char* const rhbk = "rhbk@MetaData";
  static constexpr const char* const FlagH = "FlagH@MetaData";
  static constexpr const char* const Indx = "Indx@MetaData";

  // GeoVaLs

  static constexpr const char* const geovals_orog = "surface_altitude";
  static constexpr const char* const geovals_pressure = "air_pressure";
  static constexpr const char* const geovals_pressure_rho = "air_pressure_levels";
  static constexpr const char* const geovals_height = "height";
  static constexpr const char* const geovals_height_rho = "height_levels";
  static constexpr const char* const geovals_potential_temperature = "theta";
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
    "air_temperature@DerivedObsValue";
  static constexpr const char* const eastward_wind_derived =
    "eastward_wind@DerivedObsValue";
  static constexpr const char* const northward_wind_derived =
    "northward_wind@DerivedObsValue";
  static constexpr const char* const relative_humidity_derived =
    "relative_humidity@DerivedObsValue";

  // Derived observation values (used in averaging)

  static constexpr const char* const LogP_derived = "logP@DerivedMetaData";
  static constexpr const char* const bigPgaps_derived = "bigPgaps@DerivedMetaData";

  // Derived model values (used in averaging)

  static constexpr const char* const modellevels_logP_derived = "LogPB@DerivedModelValue";
  static constexpr const char* const modellevels_ExnerP_derived = "ExnerPB@DerivedModelValue";
  static constexpr const char* const modellevels_air_temperature_derived =
    "air_temperature@DerivedModelValue";

  // Derived model values on rho levels (used in averaging)

  static constexpr const char* const modellevels_logP_rho_derived =
    "LogPA@DerivedModelValue";
  static constexpr const char* const modellevels_logPWB_rho_derived =
    "LogPWB@DerivedModelValue";
  static constexpr const char* const modellevels_ExnerP_rho_derived =
    "ExnerPA@DerivedModelValue";
};

}  // namespace ufo

#endif  // UFO_PROFILE_VARIABLENAMES_H_
