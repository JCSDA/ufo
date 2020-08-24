/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/VariableNames.h"

// Variable name: station ID

constexpr const char* const ufo::VariableNames::name_station_ID;

// Variable names: ObsValue and ObsError

constexpr const char* const ufo::VariableNames::name_obs_air_temperature;
constexpr const char* const ufo::VariableNames::name_obs_geopotential_height;
constexpr const char* const ufo::VariableNames::name_obs_eastward_wind;
constexpr const char* const ufo::VariableNames::name_obs_northward_wind;
constexpr const char* const ufo::VariableNames::name_obs_relative_humidity;
constexpr const char* const ufo::VariableNames::name_obs_dew_point_temperature;

// Variable names: HofX

constexpr const char* const ufo::VariableNames::name_hofx_air_temperature;
constexpr const char* const ufo::VariableNames::name_hofx_geopotential_height;
constexpr const char* const ufo::VariableNames::name_hofx_relative_humidity;

// Variable names: MetaData

constexpr const char* const ufo::VariableNames::name_air_pressure;
constexpr const char* const ufo::VariableNames::name_PstarBackgr;

// Variable names: QC flags

constexpr const char* const ufo::VariableNames::name_qc_ReportFlags;
constexpr const char* const ufo::VariableNames::name_qc_tFlags;
constexpr const char* const ufo::VariableNames::name_qc_RHFlags;
constexpr const char* const ufo::VariableNames::name_qc_zFlags;
constexpr const char* const ufo::VariableNames::name_qc_uFlags;
constexpr const char* const ufo::VariableNames::name_qc_vFlags;

// Variable names: counters

constexpr const char* const ufo::VariableNames::name_counter_NumAnyErrors;
constexpr const char* const ufo::VariableNames::name_counter_NumSamePErrObs;
constexpr const char* const ufo::VariableNames::name_counter_NumSuperadiabat;
constexpr const char* const ufo::VariableNames::name_counter_Num925Miss;
constexpr const char* const ufo::VariableNames::name_counter_Num100Miss;
constexpr const char* const ufo::VariableNames::name_counter_NumStdMiss;
constexpr const char* const ufo::VariableNames::name_counter_NumHydErrObs;
constexpr const char* const ufo::VariableNames::name_counter_NumIntHydErrors;
constexpr const char* const ufo::VariableNames::name_counter_NumInterpErrors;
constexpr const char* const ufo::VariableNames::name_counter_NumInterpErrObs;
constexpr const char* const ufo::VariableNames::name_counter_NumSignChange;
constexpr const char* const ufo::VariableNames::name_counter_TotCProfs;
constexpr const char* const ufo::VariableNames::name_counter_TotHProfs;
constexpr const char* const ufo::VariableNames::name_counter_TotCFlags;
constexpr const char* const ufo::VariableNames::name_counter_TotHFlags;
constexpr const char* const ufo::VariableNames::name_counter_TotLFlags;

// Variable names: corrections

constexpr const char* const ufo::VariableNames::name_tObsCorrection;
constexpr const char* const ufo::VariableNames::name_zObsCorrection;

// Variable names: intermediate values

constexpr const char* const ufo::VariableNames::name_DC;
constexpr const char* const ufo::VariableNames::name_ETol;
constexpr const char* const ufo::VariableNames::name_D;
constexpr const char* const ufo::VariableNames::name_E;
constexpr const char* const ufo::VariableNames::name_HydError;
constexpr const char* const ufo::VariableNames::name_PBottom;
constexpr const char* const ufo::VariableNames::name_StdLev;
constexpr const char* const ufo::VariableNames::name_SigAbove;
constexpr const char* const ufo::VariableNames::name_SigBelow;
constexpr const char* const ufo::VariableNames::name_IndStd;
constexpr const char* const ufo::VariableNames::name_LevErrors;
constexpr const char* const ufo::VariableNames::name_tInterp;
constexpr const char* const ufo::VariableNames::name_uInterp;
constexpr const char* const ufo::VariableNames::name_vInterp;
constexpr const char* const ufo::VariableNames::name_LogP;
constexpr const char* const ufo::VariableNames::name_NumStd;
constexpr const char* const ufo::VariableNames::name_NumSig;
constexpr const char* const ufo::VariableNames::name_Press;
constexpr const char* const ufo::VariableNames::name_Temp;
constexpr const char* const ufo::VariableNames::name_rh;
constexpr const char* const ufo::VariableNames::name_td;
constexpr const char* const ufo::VariableNames::name_tbk;
constexpr const char* const ufo::VariableNames::name_rhbk;
constexpr const char* const ufo::VariableNames::name_FlagH;
constexpr const char* const ufo::VariableNames::name_Indx;
