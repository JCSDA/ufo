/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/VariableNames.h"

// Variable names: observation values

constexpr const char* const ufo::VariableNames::obs_air_pressure;
constexpr const char* const ufo::VariableNames::obs_air_temperature;
constexpr const char* const ufo::VariableNames::obs_relative_humidity;
constexpr const char* const ufo::VariableNames::obs_eastward_wind;
constexpr const char* const ufo::VariableNames::obs_northward_wind;
constexpr const char* const ufo::VariableNames::obs_geopotential_height;
constexpr const char* const ufo::VariableNames::obs_dew_point_temperature;

// Variable names: observation errors

constexpr const char* const ufo::VariableNames::obserr_air_temperature;
constexpr const char* const ufo::VariableNames::obserr_relative_humidity;
constexpr const char* const ufo::VariableNames::obserr_eastward_wind;
constexpr const char* const ufo::VariableNames::obserr_northward_wind;
constexpr const char* const ufo::VariableNames::obserr_geopotential_height;
constexpr const char* const ufo::VariableNames::obserr_dew_point_temperature;

// Variable names: HofX

constexpr const char* const ufo::VariableNames::hofx_air_temperature;
constexpr const char* const ufo::VariableNames::hofx_geopotential_height;
constexpr const char* const ufo::VariableNames::hofx_relative_humidity;
constexpr const char* const ufo::VariableNames::hofx_eastward_wind;
constexpr const char* const ufo::VariableNames::hofx_northward_wind;
constexpr const char* const ufo::VariableNames::hofx_dew_point_temperature;

// Variable names: background errors

constexpr const char* const ufo::VariableNames::bkgerr_air_temperature;
constexpr const char* const ufo::VariableNames::bkgerr_relative_humidity;
constexpr const char* const ufo::VariableNames::bkgerr_eastward_wind;
constexpr const char* const ufo::VariableNames::bkgerr_northward_wind;
constexpr const char* const ufo::VariableNames::bkgerr_geopotential_height;
constexpr const char* const ufo::VariableNames::bkgerr_dew_point_temperature;

// Variable names: probability of gross error

constexpr const char* const ufo::VariableNames::pge_air_temperature;
constexpr const char* const ufo::VariableNames::pge_relative_humidity;
constexpr const char* const ufo::VariableNames::pge_eastward_wind;
constexpr const char* const ufo::VariableNames::pge_northward_wind;
constexpr const char* const ufo::VariableNames::pge_geopotential_height;

// Variable names: probability of gross error used in buddy check

constexpr const char* const ufo::VariableNames::pgebd_air_temperature;
constexpr const char* const ufo::VariableNames::pgebd_relative_humidity;
constexpr const char* const ufo::VariableNames::pgebd_eastward_wind;
constexpr const char* const ufo::VariableNames::pgebd_northward_wind;
constexpr const char* const ufo::VariableNames::pgebd_geopotential_height;

// Variable names: MetaData

constexpr const char* const ufo::VariableNames::station_ID;
constexpr const char* const ufo::VariableNames::obs_level_time;
constexpr const char* const ufo::VariableNames::PstarBackgr;
constexpr const char* const ufo::VariableNames::ObsType;
constexpr const char* const ufo::VariableNames::Latitude;
constexpr const char* const ufo::VariableNames::Longitude;
constexpr const char* const ufo::VariableNames::Time;
constexpr const char* const ufo::VariableNames::Zstation;

// Variable names: QC flags

constexpr const char* const ufo::VariableNames::qcflags_observation_report;
constexpr const char* const ufo::VariableNames::qcflags_air_temperature;
constexpr const char* const ufo::VariableNames::qcflags_relative_humidity;
constexpr const char* const ufo::VariableNames::qcflags_geopotential_height;
constexpr const char* const ufo::VariableNames::qcflags_eastward_wind;
constexpr const char* const ufo::VariableNames::qcflags_northward_wind;
constexpr const char* const ufo::VariableNames::qcflags_time;

// Variable names: counters

constexpr const char* const ufo::VariableNames::counter_NumAnyErrors;
constexpr const char* const ufo::VariableNames::counter_NumSamePErrObs;
constexpr const char* const ufo::VariableNames::counter_NumSuperadiabat;
constexpr const char* const ufo::VariableNames::counter_Num925Miss;
constexpr const char* const ufo::VariableNames::counter_Num100Miss;
constexpr const char* const ufo::VariableNames::counter_NumStdMiss;
constexpr const char* const ufo::VariableNames::counter_NumHydErrObs;
constexpr const char* const ufo::VariableNames::counter_NumIntHydErrors;
constexpr const char* const ufo::VariableNames::counter_NumInterpErrors;
constexpr const char* const ufo::VariableNames::counter_NumInterpErrObs;
constexpr const char* const ufo::VariableNames::counter_NumSignChange;
constexpr const char* const ufo::VariableNames::counter_TotCProfs;
constexpr const char* const ufo::VariableNames::counter_TotHProfs;
constexpr const char* const ufo::VariableNames::counter_TotCFlags;
constexpr const char* const ufo::VariableNames::counter_TotHFlags;
constexpr const char* const ufo::VariableNames::counter_TotLFlags;

// Variable names: corrections

constexpr const char* const ufo::VariableNames::obscorrection_air_temperature;
constexpr const char* const ufo::VariableNames::obscorrection_geopotential_height;

// Variable names: intermediate values

constexpr const char* const ufo::VariableNames::DC;
constexpr const char* const ufo::VariableNames::ETol;
constexpr const char* const ufo::VariableNames::D;
constexpr const char* const ufo::VariableNames::E;
constexpr const char* const ufo::VariableNames::HydError;
constexpr const char* const ufo::VariableNames::PBottom;
constexpr const char* const ufo::VariableNames::StdLev;
constexpr const char* const ufo::VariableNames::SigAbove;
constexpr const char* const ufo::VariableNames::SigBelow;
constexpr const char* const ufo::VariableNames::IndStd;
constexpr const char* const ufo::VariableNames::LevErrors;
constexpr const char* const ufo::VariableNames::tInterp;
constexpr const char* const ufo::VariableNames::uInterp;
constexpr const char* const ufo::VariableNames::vInterp;
constexpr const char* const ufo::VariableNames::LogP;
constexpr const char* const ufo::VariableNames::NumStd;
constexpr const char* const ufo::VariableNames::NumSig;
constexpr const char* const ufo::VariableNames::Press;
constexpr const char* const ufo::VariableNames::Temp;
constexpr const char* const ufo::VariableNames::rh;
constexpr const char* const ufo::VariableNames::td;
constexpr const char* const ufo::VariableNames::tbk;
constexpr const char* const ufo::VariableNames::rhbk;
constexpr const char* const ufo::VariableNames::FlagH;
constexpr const char* const ufo::VariableNames::Indx;
