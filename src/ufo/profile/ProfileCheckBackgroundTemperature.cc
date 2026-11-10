/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckBackgroundTemperature.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundTemperature>
  makerProfileCheckBackgroundTemperature_("BackgroundTemperature");

  ProfileCheckBackgroundTemperature::ProfileCheckBackgroundTemperature
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBackgroundTemperature::doPrepareOutputColumns(
      ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Background check for temperature" << std::endl;

    const std::vector<std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::pge_air_temperature,
       ufo::ProfileVariableNames::obscorrection_air_temperature};
    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_superadiabat_t,
       ufo::ProfileVariableNames::diagflags_interpolation_t,
       ufo::ProfileVariableNames::diagflags_hydro_t,
       ufo::ProfileVariableNames::diagflags_back_perf_t,
       ufo::ProfileVariableNames::diagflags_back_reject_t,
       ufo::ProfileVariableNames::diagflags_perm_reject_t,
       ufo::ProfileVariableNames::diagflags_final_reject_t};

    for (const auto &variableFloat : variableNamesFloat) {
      profileDataHandler.get<float>(variableFloat);
    }
    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileCheckBackgroundTemperature::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <float> &tPGE =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::pge_air_temperature);
    std::vector <bool> &diagFlagsTSuperadiabat =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_superadiabat_t);
    std::vector <bool> &diagFlagsTInterp =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_t);
    std::vector <bool> &diagFlagsTHydro =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_hydro_t);
    std::vector <bool> &diagFlagsTBackPerf =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_perf_t);
    std::vector <bool> &diagFlagsTBackReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_reject_t);
    std::vector <bool> &diagFlagsTPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_t);
    std::vector <bool> &diagFlagsTFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    const std::vector <float> &tObsCorrection =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obscorrection_air_temperature);
    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &Latitude =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::Latitude);
    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_temperature);
    const std::vector <float> &tObsErr =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obserr_air_temperature);
    const std::vector <float> &tBkg =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_air_temperature);
    const std::vector <float> &tBkgErr =
      profileDataHandler.get<float>
      (options_.bkgErrGroup.value() + "/" + options_.bkgErrName_air_temperature.value());
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(Latitude, pressures,
                                         tObs, tObsErr, tBkg, tBkgErr,
                                         tPGE,
                                         diagFlagsTSuperadiabat,
                                         diagFlagsTInterp,
                                         diagFlagsTHydro,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(Latitude, pressures,
                                                    tObs, tObsErr, tBkg, tBkgErr,
                                                    tPGE,
                                                    diagFlagsTSuperadiabat,
                                                    diagFlagsTInterp,
                                                    diagFlagsTHydro,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, options_.BkCheck_PdBad_t.value());
    // Local version of temperature background error.
    std::vector <float> BackgrErrT(numProfileLevels, 0);
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      BackgrErrT[jlev] = tBkgErr[jlev];
    }
    // Extra representivity error for data on reported levels.
    if (!ModelLevels) {
      const float Psplit =
        std::fabs(Latitude[0]) < options_.BkCheck_Psplit_latitude_tropics ?
                                 options_.BkCheck_Psplit_tropics.value() :
                                 options_.BkCheck_Psplit_extratropics.value();
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (pressures[jlev] <= Psplit) {
          BackgrErrT[jlev] = tBkgErr[jlev] == missingValueFloat ?
            missingValueFloat :
            tBkgErr[jlev] * options_.BkCheck_ErrorInflationBelowPsplit.value();
        } else {
          BackgrErrT[jlev] = tBkgErr[jlev] == missingValueFloat ?
            missingValueFloat :
            tBkgErr[jlev] * options_.BkCheck_ErrorInflationAbovePsplit.value();
        }
      }
    }

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (diagFlagsTSuperadiabat[jlev])
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
      if (diagFlagsTInterp[jlev])
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
      if (diagFlagsTHydro[jlev])
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           tObsFinal,
                           tObsErr,
                           tBkg,
                           BackgrErrT,  // Used instead of tBkgErr.
                           PdBad,
                           ModelLevels,
                           diagFlagsTBackPerf,
                           diagFlagsTBackReject,
                           diagFlagsTPermReject,
                           diagFlagsTFinalReject,
                           tPGE);
  }
}  // namespace ufo
