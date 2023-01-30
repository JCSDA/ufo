/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/VariableNames.h"

// Observation values

constexpr const char* const ufo::VariableNames::obs_air_pressure;
constexpr const char* const ufo::VariableNames::obs_air_temperature;
constexpr const char* const ufo::VariableNames::obs_relative_humidity;
constexpr const char* const ufo::VariableNames::obs_eastward_wind;
constexpr const char* const ufo::VariableNames::obs_northward_wind;
constexpr const char* const ufo::VariableNames::obs_geopotential_height;
constexpr const char* const ufo::VariableNames::obs_dew_point_temperature;

// Derived observation values

constexpr const char* const ufo::VariableNames::obs_derived_air_pressure;

// Observation errors

constexpr const char* const ufo::VariableNames::obserr_air_temperature;
constexpr const char* const ufo::VariableNames::obserr_relative_humidity;
constexpr const char* const ufo::VariableNames::obserr_eastward_wind;
constexpr const char* const ufo::VariableNames::obserr_northward_wind;
constexpr const char* const ufo::VariableNames::obserr_geopotential_height;
constexpr const char* const ufo::VariableNames::obserr_dew_point_temperature;

// HofX

constexpr const char* const ufo::VariableNames::hofx_air_temperature;
constexpr const char* const ufo::VariableNames::hofx_geopotential_height;
constexpr const char* const ufo::VariableNames::hofx_relative_humidity;
constexpr const char* const ufo::VariableNames::hofx_eastward_wind;
constexpr const char* const ufo::VariableNames::hofx_northward_wind;
constexpr const char* const ufo::VariableNames::hofx_dew_point_temperature;

// Probability of gross error

constexpr const char* const ufo::VariableNames::pge_air_temperature;
constexpr const char* const ufo::VariableNames::pge_relative_humidity;
constexpr const char* const ufo::VariableNames::pge_eastward_wind;
constexpr const char* const ufo::VariableNames::pge_northward_wind;
constexpr const char* const ufo::VariableNames::pge_geopotential_height;

// MetaData

constexpr const char* const ufo::VariableNames::station_ID;
constexpr const char* const ufo::VariableNames::ObsType;
constexpr const char* const ufo::VariableNames::Latitude;
constexpr const char* const ufo::VariableNames::Longitude;
constexpr const char* const ufo::VariableNames::Zstation;
constexpr const char* const ufo::VariableNames::LevelType;
constexpr const char* const ufo::VariableNames::InstrType;
constexpr const char* const ufo::VariableNames::extended_obs_space;

// QC flags

constexpr const char* const ufo::VariableNames::qcflags_observation_report;
constexpr const char* const ufo::VariableNames::qcflags_air_temperature;
constexpr const char* const ufo::VariableNames::qcflags_relative_humidity;
constexpr const char* const ufo::VariableNames::qcflags_geopotential_height;
constexpr const char* const ufo::VariableNames::qcflags_eastward_wind;
constexpr const char* const ufo::VariableNames::qcflags_northward_wind;
constexpr const char* const ufo::VariableNames::qcflags_wind_profiler;

// Diagnostic flags

constexpr const char* const ufo::VariableNames::diagflags_profile_interpolation_eastward_wind;
constexpr const char* const ufo::VariableNames::diagflags_profile_standard_level_eastward_wind;

// Counters

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
constexpr const char* const ufo::VariableNames::counter_NumGapsT;
constexpr const char* const ufo::VariableNames::counter_NumGapsU;
constexpr const char* const ufo::VariableNames::counter_NumGapsUWP;
constexpr const char* const ufo::VariableNames::counter_NumGapsRH;

// Corrections

constexpr const char* const ufo::VariableNames::obscorrection_air_temperature;
constexpr const char* const ufo::VariableNames::obscorrection_geopotential_height;

// Intermediate values

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

// GeoVaLs

constexpr const char* const ufo::VariableNames::geovals_orog;
constexpr const char* const ufo::VariableNames::geovals_pressure;
constexpr const char* const ufo::VariableNames::geovals_pressure_rho;
constexpr const char* const ufo::VariableNames::geovals_pressure_rho_minus_one;
constexpr const char* const ufo::VariableNames::geovals_height;
constexpr const char* const ufo::VariableNames::geovals_height_rho;
constexpr const char* const ufo::VariableNames::geovals_height_rho_minus_one;
constexpr const char* const ufo::VariableNames::geovals_potential_temperature;
constexpr const char* const ufo::VariableNames::geovals_air_temperature;
constexpr const char* const ufo::VariableNames::geovals_surface_pressure;
constexpr const char* const ufo::VariableNames::geovals_relative_humidity;

// GeoVaLs used in validation

constexpr const char* const ufo::VariableNames::geovals_testreference_logP;
constexpr const char* const ufo::VariableNames::geovals_testreference_ExnerP;
constexpr const char* const ufo::VariableNames::geovals_testreference_logP_rho;
constexpr const char* const ufo::VariableNames::geovals_testreference_ExnerP_rho;
constexpr const char* const ufo::VariableNames::geovals_testreference_air_temperature;
constexpr const char* const ufo::VariableNames::geovals_testreference_eastward_wind;
constexpr const char* const ufo::VariableNames::geovals_testreference_northward_wind;
constexpr const char* const ufo::VariableNames::geovals_testreference_relative_humidity;
constexpr const char* const ufo::VariableNames::geovals_testreference_air_temperature_qcflags;
constexpr const char* const ufo::VariableNames::geovals_testreference_eastward_wind_qcflags;
constexpr const char* const ufo::VariableNames::geovals_testreference_northward_wind_qcflags;
constexpr const char* const ufo::VariableNames::geovals_testreference_relative_humidity_qcflags;

// Averaged values on model levels

constexpr const char* const ufo::VariableNames::air_temperature_derived;
constexpr const char* const ufo::VariableNames::eastward_wind_derived;
constexpr const char* const ufo::VariableNames::northward_wind_derived;
constexpr const char* const ufo::VariableNames::relative_humidity_derived;

// Derived observation values (used in averaging)

constexpr const char* const ufo::VariableNames::LogP_derived;
constexpr const char* const ufo::VariableNames::bigPgaps_derived;

// Derived model values (used in averaging)

constexpr const char* const ufo::VariableNames::modellevels_logP_derived;
constexpr const char* const ufo::VariableNames::modellevels_ExnerP_derived;
constexpr const char* const ufo::VariableNames::modellevels_air_temperature_derived;

// Derived model values on rho levels (used in averaging)

constexpr const char* const ufo::VariableNames::modellevels_logP_rho_derived;
constexpr const char* const ufo::VariableNames::modellevels_logPWB_rho_derived;
constexpr const char* const ufo::VariableNames::modellevels_ExnerP_rho_derived;
