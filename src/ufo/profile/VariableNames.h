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
  static constexpr const char* const obs_geopotential_height = "geopotential_height@ObsValue";
  static constexpr const char* const obs_dew_point_temperature = "dew_point_temperature@ObsValue";

  // Observation errors

  static constexpr const char* const obserr_air_temperature = "air_temperature@ObsError";
  static constexpr const char* const obserr_relative_humidity = "relative_humidity@ObsError";
  static constexpr const char* const obserr_eastward_wind = "eastward_wind@ObsError";
  static constexpr const char* const obserr_northward_wind = "northward_wind@ObsError";
  static constexpr const char* const obserr_geopotential_height = "geopotential_height@ObsError";
  static constexpr const char* const obserr_dew_point_temperature =
    "dew_point_temperature@ObsError";

  // HofX

  static constexpr const char* const hofx_air_temperature = "air_temperature@HofX";
  static constexpr const char* const hofx_geopotential_height = "geopotential_height@HofX";
  static constexpr const char* const hofx_relative_humidity = "relative_humidity@HofX";
  static constexpr const char* const hofx_eastward_wind = "eastward_wind@HofX";
  static constexpr const char* const hofx_northward_wind = "northward_wind@HofX";
  static constexpr const char* const hofx_dew_point_temperature = "dew_point_temperature@HofX";

  // Background errors

  static constexpr const char* const bkgerr_air_temperature = "air_temperature@BkgError";
  static constexpr const char* const bkgerr_relative_humidity = "relative_humidity@BkgError";
  static constexpr const char* const bkgerr_eastward_wind = "eastward_wind@BkgError";
  static constexpr const char* const bkgerr_northward_wind = "northward_wind@BkgError";
  static constexpr const char* const bkgerr_geopotential_height = "geopotential_height@BkgError";
  static constexpr const char* const bkgerr_dew_point_temperature =
    "dew_point_temperature@BkgError";

  // Probability of gross error

  static constexpr const char* const pge_air_temperature = "air_temperature@GrossErrorProbability";
  static constexpr const char* const pge_relative_humidity =
    "relative_humidity@GrossErrorProbability";
  static constexpr const char* const pge_eastward_wind = "eastward_wind@GrossErrorProbability";
  static constexpr const char* const pge_northward_wind = "northward_wind@GrossErrorProbability";
  static constexpr const char* const pge_geopotential_height =
    "geopotential_height@GrossErrorProbability";

  // Probability of gross error used in buddy check

  static constexpr const char* const pgebd_air_temperature =
    "air_temperature@GrossErrorProbabilityBuddyCheck";
  static constexpr const char* const pgebd_relative_humidity =
    "relative_humidity@GrossErrorProbabilityBuddyCheck";
  static constexpr const char* const pgebd_eastward_wind =
    "eastward_wind@GrossErrorProbabilityBuddyCheck";
  static constexpr const char* const pgebd_northward_wind =
    "northward_wind@GrossErrorProbabilityBuddyCheck";
  static constexpr const char* const pgebd_geopotential_height =
    "geopotential_height@GrossErrorProbabilityBuddyCheck";

  // MetaData

  static constexpr const char* const station_ID = "station_id@MetaData";
  static constexpr const char* const obs_level_time = "level_time@MetaData";
  static constexpr const char* const ObsType = "ObsType@MetaData";
  static constexpr const char* const Latitude = "latitude@MetaData";
  static constexpr const char* const Longitude = "longitude@MetaData";
  static constexpr const char* const Time = "time@MetaData";
  static constexpr const char* const Zstation = "Zstation@MetaData";
  static constexpr const char* const LevelType = "LevelType@MetaData";
  static constexpr const char* const InstrType = "InstrType@MetaData";

  // QC flags

  static constexpr const char* const qcflags_observation_report = "observation_report@QCFlags";
  static constexpr const char* const qcflags_air_temperature = "air_temperature@QCFlags";
  static constexpr const char* const qcflags_relative_humidity = "relative_humidity@QCFlags";
  static constexpr const char* const qcflags_geopotential_height = "geopotential_height@QCFlags";
  static constexpr const char* const qcflags_eastward_wind = "eastward_wind@QCFlags";
  static constexpr const char* const qcflags_northward_wind = "northward_wind@QCFlags";
  static constexpr const char* const qcflags_time = "time@QCFlags";
  static constexpr const char* const qcflags_wind_profiler = "wind_profiler@QCFlags";

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
    "geopotential_height@Corrections";

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

  // Derived values

  static constexpr const char* const LogP_derived = "logP@DerivedValue";
  static constexpr const char* const bigPgaps_derived = "bigPgaps@DerivedValue";

  // GeoVaLs

  static constexpr const char* const geovals_orog = "orography";
  static constexpr const char* const geovals_pressure = "air_pressure";
  static constexpr const char* const geovals_pressure_rho = "air_pressure_rho";
  static constexpr const char* const geovals_height = "height";
  static constexpr const char* const geovals_potential_temperature = "air_potential_temperature";
  static constexpr const char* const geovals_surface_pressure = "surface_pressure";
  static constexpr const char* const geovals_relative_humidity = "relative_humidity";

  // GeoVaLs used in validation

  static constexpr const char* const geovals_logP = "LogPB";
  static constexpr const char* const geovals_ExnerP = "ExnerPB";
  static constexpr const char* const geovals_logP_rho = "LogPA";
  static constexpr const char* const geovals_ExnerP_rho = "ExnerPA";
  static constexpr const char* const geovals_air_temperature = "air_temperature";
  static constexpr const char* const geovals_average_air_temperature = "average_air_temperature";
  static constexpr const char* const geovals_average_eastward_wind = "average_eastward_wind";
  static constexpr const char* const geovals_average_northward_wind = "average_northward_wind";
  static constexpr const char* const geovals_average_relative_humidity =
    "average_relative_humidity";
  static constexpr const char* const geovals_average_air_temperature_qcflags =
    "average_air_temperature_flags";
  static constexpr const char* const geovals_average_eastward_wind_qcflags =
    "average_eastward_wind_flags";
  static constexpr const char* const geovals_average_northward_wind_qcflags =
    "average_northward_wind_flags";
  static constexpr const char* const geovals_average_relative_humidity_qcflags =
    "average_relative_humidity_flags";

  // Derived values on model levels

  static constexpr const char* const geovals_logP_derived = "LogPB@ModelLevelsDerivedValue";
  static constexpr const char* const geovals_ExnerP_derived = "ExnerPB@ModelLevelsDerivedValue";
  static constexpr const char* const geovals_air_temperature_derived =
    "air_temperature@ModelLevelsDerivedValue";

  // Averaged values on model levels

  static constexpr const char* const average_air_temperature_derived =
    "average_air_temperature@ModelLevelsDerivedValue";
  static constexpr const char* const average_eastward_wind_derived =
    "average_eastward_wind@ModelLevelsDerivedValue";
  static constexpr const char* const average_northward_wind_derived =
    "average_northward_wind@ModelLevelsDerivedValue";
  static constexpr const char* const average_relative_humidity_derived =
    "average_relative_humidity@ModelLevelsDerivedValue";
  static constexpr const char* const average_air_temperature_qcflags =
    "average_air_temperature@ModelLevelsQCFlags";
  static constexpr const char* const average_eastward_wind_qcflags =
    "average_eastward_wind@ModelLevelsQCFlags";
  static constexpr const char* const average_northward_wind_qcflags =
    "average_northward_wind@ModelLevelsQCFlags";
  static constexpr const char* const average_relative_humidity_qcflags =
    "average_relative_humidity@ModelLevelsQCFlags";

  // Derived values on model rho levels

  static constexpr const char* const geovals_logP_rho_derived = "LogPA@ModelRhoLevelsDerivedValue";
  static constexpr const char* const geovals_logPWB_rho_derived =
    "LogPWB@ModelRhoLevelsDerivedValue";
  static constexpr const char* const geovals_ExnerP_rho_derived =
    "ExnerPA@ModelRhoLevelsDerivedValue";
};

}  // namespace ufo

#endif  // UFO_PROFILE_VARIABLENAMES_H_
