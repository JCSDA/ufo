/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckBasic.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBasic> makerProfileCheckBasic_("Basic");

  ProfileCheckBasic::ProfileCheckBasic(const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBasic::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    if (options_.BChecks_Skip.value())
      return;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_t,
       ufo::ProfileVariableNames::diagflags_final_reject_u,
       ufo::ProfileVariableNames::diagflags_final_reject_v,
       ufo::ProfileVariableNames::diagflags_final_reject_rh};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileCheckBasic::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Basic checks" << std::endl;

    // Set basic check result to true
    result_ = true;

    // Skip this routine if specifically requested
    if (options_.BChecks_Skip.value())
      {
        oops::Log::debug() << "Skipping basic checks" << std::endl;
        return;
      }

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    // All diagnostic flags are retrieved for the basic checks.
    // (Some might be empty; that is checked before they are used.)
    std::vector <bool> &diagFlagsTFinalReject = profileDataHandler.get<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <bool> &diagFlagsUFinalReject = profileDataHandler.get<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_u);
    std::vector <bool> &diagFlagsVFinalReject = profileDataHandler.get<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_v);
    std::vector <bool> &diagFlagsRHFinalReject = profileDataHandler.get<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_rh);

    // Warn and exit if pressures vector is empty
    if (pressures.empty()) {
      result_ = false;
      oops::Log::debug() << "Pressures vector is empty" << std::endl;
      return;
     }

    // Apply a series of basic checks and print the results to the debug stream.

    // Is the number of levels to check OK?
    bool numProfileLevelsOK = (numProfileLevels > 0);

    // Are any levels in the wrong order?
    bool pressOrderOK = true;
    for (int jlev = 0; jlev < numProfileLevels - 1; ++jlev) {
      pressOrderOK = pressOrderOK && (pressures[jlev] >= pressures[jlev + 1]);
      if (!pressOrderOK) break;
    }

    // Is the pressure at the first level > maximum value pressure?
    bool maxPressOK = (pressures.size() > 0 ?
                       pressures.front() <= options_.BChecks_maxValidP.value() :
                       false);

    // Is the pressure at the final level < minimum value pressure?
    bool minPressOK = (pressures.size() > 0 ?
                       pressures.back() > options_.BChecks_minValidP.value() :
                       false);

    oops::Log::debug() << " -> numProfileLevelsOK: " << numProfileLevelsOK << std::endl;
    oops::Log::debug() << " -> pressOrderOK: " << pressOrderOK << std::endl;
    oops::Log::debug() << " -> maxPressOK: " << maxPressOK << std::endl;
    oops::Log::debug() << " -> minPressOK: " << minPressOK << std::endl;

    // If `result_` is true then no further profile QC checks instantiated by this filter will run.
    result_ = numProfileLevelsOK && pressOrderOK && maxPressOK && minPressOK;
    oops::Log::debug() << " -> basicResult: " << result_ << std::endl;

    // Flag any failing levels if required.
    if (options_.flagBasicChecksFail) {
      // Starting value of pressure (Pa).
      float pressure_current = options_.BChecks_maxValidP.value();
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (pressures[jlev] == missingValueFloat) continue;
        if (pressures[jlev] > pressure_current ||
            pressures[jlev] <= options_.BChecks_minValidP.value()) {
          if (!diagFlagsTFinalReject.empty()) diagFlagsTFinalReject[jlev] = true;
          if (!diagFlagsUFinalReject.empty()) diagFlagsUFinalReject[jlev] = true;
          if (!diagFlagsVFinalReject.empty()) diagFlagsVFinalReject[jlev] = true;
          if (!diagFlagsRHFinalReject.empty()) diagFlagsRHFinalReject[jlev] = true;
        }
        pressure_current = pressures[jlev];
      }
    }
  }
}  // namespace ufo

