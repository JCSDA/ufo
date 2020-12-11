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
  // Variable names: observation values

  static constexpr const char* const obs_air_pressure = "air_pressure@MetaData";
  static constexpr const char* const obs_air_temperature = "air_temperature@ObsValue";
  static constexpr const char* const obs_relative_humidity = "relative_humidity@ObsValue";
  static constexpr const char* const obs_eastward_wind = "eastward_wind@ObsValue";
  static constexpr const char* const obs_northward_wind = "northward_wind@ObsValue";
  static constexpr const char* const obs_geopotential_height = "geopotential_height@ObsValue";
  static constexpr const char* const obs_dew_point_temperature = "dew_point_temperature@ObsValue";

  // Variable names: observation errors

  static constexpr const char* const obserr_air_temperature = "air_temperature@ObsError";
  static constexpr const char* const obserr_relative_humidity = "relative_humidity@ObsError";
  static constexpr const char* const obserr_eastward_wind = "eastward_wind@ObsError";
  static constexpr const char* const obserr_northward_wind = "northward_wind@ObsError";
  static constexpr const char* const obserr_geopotential_height = "geopotential_height@ObsError";
  static constexpr const char* const obserr_dew_point_temperature =
    "dew_point_temperature@ObsError";

  // Variable names: HofX

  static constexpr const char* const hofx_air_temperature = "air_temperature@HofX";
  static constexpr const char* const hofx_geopotential_height = "geopotential_height@HofX";
  static constexpr const char* const hofx_relative_humidity = "relative_humidity@HofX";
  static constexpr const char* const hofx_eastward_wind = "eastward_wind@HofX";
  static constexpr const char* const hofx_northward_wind = "northward_wind@HofX";
  static constexpr const char* const hofx_dew_point_temperature = "dew_point_temperature@HofX";

  // Variable names: background errors

  static constexpr const char* const bkgerr_air_temperature = "air_temperature@BkgError";
  static constexpr const char* const bkgerr_relative_humidity = "relative_humidity@BkgError";
  static constexpr const char* const bkgerr_eastward_wind = "eastward_wind@BkgError";
  static constexpr const char* const bkgerr_northward_wind = "northward_wind@BkgError";
  static constexpr const char* const bkgerr_geopotential_height = "geopotential_height@BkgError";
  static constexpr const char* const bkgerr_dew_point_temperature =
    "dew_point_temperature@BkgError";

  // Variable names: probability of gross error

  static constexpr const char* const pge_air_temperature = "air_temperature@GrossErrorProbability";
  static constexpr const char* const pge_relative_humidity =
    "relative_humidity@GrossErrorProbability";
  static constexpr const char* const pge_eastward_wind = "eastward_wind@GrossErrorProbability";
  static constexpr const char* const pge_northward_wind = "northward_wind@GrossErrorProbability";
  static constexpr const char* const pge_geopotential_height =
    "geopotential_height@GrossErrorProbability";

  // Variable names: probability of gross error used in buddy check

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

  // Variable names: MetaData

  static constexpr const char* const station_ID = "station_id@MetaData";
  static constexpr const char* const obs_level_time = "level_time@MetaData";
  static constexpr const char* const PstarBackgr = "PstarBackgr@MetaData";
  static constexpr const char* const ObsType = "ObsType@MetaData";
  static constexpr const char* const Latitude = "latitude@MetaData";
  static constexpr const char* const Longitude = "longitude@MetaData";
  static constexpr const char* const Time = "time@MetaData";
  static constexpr const char* const Zstation = "Zstation@MetaData";
  static constexpr const char* const LevelType = "LevelType@MetaData";

  // Variable names: QC flags

  static constexpr const char* const qcflags_observation_report = "observation_report@QCFlags";
  static constexpr const char* const qcflags_air_temperature = "air_temperature@QCFlags";
  static constexpr const char* const qcflags_relative_humidity = "relative_humidity@QCFlags";
  static constexpr const char* const qcflags_geopotential_height = "geopotential_height@QCFlags";
  static constexpr const char* const qcflags_eastward_wind = "eastward_wind@QCFlags";
  static constexpr const char* const qcflags_northward_wind = "northward_wind@QCFlags";
  static constexpr const char* const qcflags_time = "time@QCFlags";
  static constexpr const char* const qcflags_wind_profiler = "wind_profiler@QCFlags";

  // Variable names: counters

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

  // Variable names: corrections

  static constexpr const char* const obscorrection_air_temperature = "air_temperature@Corrections";
  static constexpr const char* const obscorrection_geopotential_height =
    "geopotential_height@Corrections";

  // Variable names: intermediate values

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
};

}  // namespace ufo

#endif  // UFO_PROFILE_VARIABLENAMES_H_
