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
  // Variable name: station ID

  static constexpr const char* const name_station_ID = "station_id@MetaData";

  // Variable names: ObsValue and ObsError

  static constexpr const char* const name_obs_air_temperature = "air_temperature@ObsValue";
  static constexpr const char* const name_obs_geopotential_height = "geopotential_height@ObsValue";
  static constexpr const char* const name_obs_eastward_wind = "eastward_wind@ObsValue";
  static constexpr const char* const name_obs_northward_wind = "northward_wind@ObsValue";
  static constexpr const char* const name_obs_relative_humidity = "relative_humidity@ObsValue";
  static constexpr const char* const name_obs_dew_point_temperature =
    "dew_point_temperature@ObsValue";

  // Variable names: HofX

  static constexpr const char* const name_hofx_air_temperature = "air_temperature@HofX";
  static constexpr const char* const name_hofx_geopotential_height = "geopotential_height@HofX";
  static constexpr const char* const name_hofx_relative_humidity = "relative_humidity@HofX";

  // Variable names: MetaData

  static constexpr const char* const name_air_pressure = "air_pressure@MetaData";
  static constexpr const char* const name_PstarBackgr = "PstarBackgr@MetaData";

  // Variable names: QC flags

  static constexpr const char* const name_qc_ReportFlags = "ReportFlags@QCFlags";
  static constexpr const char* const name_qc_tFlags = "tFlags@QCFlags";
  static constexpr const char* const name_qc_RHFlags = "RHFlags@QCFlags";
  static constexpr const char* const name_qc_zFlags = "zFlags@QCFlags";
  static constexpr const char* const name_qc_uFlags = "uFlags@QCFlags";
  static constexpr const char* const name_qc_vFlags = "vFlags@QCFlags";

  // Variable names: counters

  static constexpr const char* const name_counter_NumAnyErrors = "NumAnyErrors@Counters";
  static constexpr const char* const name_counter_NumSamePErrObs = "NumSamePErrObs@Counters";
  static constexpr const char* const name_counter_NumSuperadiabat = "NumSuperadiabat@Counters";
  static constexpr const char* const name_counter_Num925Miss = "Num925Miss@Counters";
  static constexpr const char* const name_counter_Num100Miss = "Num100Miss@Counters";
  static constexpr const char* const name_counter_NumStdMiss = "NumStdMiss@Counters";
  static constexpr const char* const name_counter_NumHydErrObs = "NumHydErrObs@Counters";
  static constexpr const char* const name_counter_NumIntHydErrors = "NumIntHydErrors@Counters";
  static constexpr const char* const name_counter_NumInterpErrors = "NumInterpErrors@Counters";
  static constexpr const char* const name_counter_NumInterpErrObs = "NumInterpErrObs@Counters";
  static constexpr const char* const name_counter_NumSignChange = "NumSignChange@Counters";
  static constexpr const char* const name_counter_TotCProfs = "TotCProfs@Counters";
  static constexpr const char* const name_counter_TotHProfs = "TotHProfs@Counters";
  static constexpr const char* const name_counter_TotCFlags = "TotCFlags@Counters";
  static constexpr const char* const name_counter_TotHFlags = "TotHFlags@Counters";
  static constexpr const char* const name_counter_TotLFlags = "TotLFlags@Counters";

  // Variable names: corrections

  static constexpr const char* const name_tObsCorrection = "tObsCorrection@Corrections";
  static constexpr const char* const name_zObsCorrection = "zObsCorrection@Corrections";

  // Variable names: intermediate values

  static constexpr const char* const name_DC = "DC@MetaData";
  static constexpr const char* const name_ETol = "ETol@MetaData";
  static constexpr const char* const name_D = "D@MetaData";
  static constexpr const char* const name_E = "E@MetaData";
  static constexpr const char* const name_HydError = "HydError@MetaData";
  static constexpr const char* const name_PBottom = "PBottom@MetaData";
  static constexpr const char* const name_StdLev = "StdLev@MetaData";
  static constexpr const char* const name_SigAbove = "SigAbove@MetaData";
  static constexpr const char* const name_SigBelow = "SigBelow@MetaData";
  static constexpr const char* const name_IndStd = "IndStd@MetaData";
  static constexpr const char* const name_LevErrors = "LevErrors@MetaData";
  static constexpr const char* const name_tInterp = "tInterp@MetaData";
  static constexpr const char* const name_uInterp = "uInterp@MetaData";
  static constexpr const char* const name_vInterp = "vInterp@MetaData";
  static constexpr const char* const name_LogP = "LogP@MetaData";
  static constexpr const char* const name_NumStd = "NumStd@MetaData";
  static constexpr const char* const name_NumSig = "NumSig@MetaData";
  static constexpr const char* const name_Press = "Press@MetaData";
  static constexpr const char* const name_Temp = "Temp@MetaData";
  static constexpr const char* const name_rh = "rh@MetaData";
  static constexpr const char* const name_td = "td@MetaData";
  static constexpr const char* const name_tbk = "tbk@MetaData";
  static constexpr const char* const name_rhbk = "rhbk@MetaData";
  static constexpr const char* const name_FlagH = "FlagH@MetaData";
  static constexpr const char* const name_Indx = "Indx@MetaData";
};

}  // namespace ufo

#endif  // UFO_PROFILE_VARIABLENAMES_H_
