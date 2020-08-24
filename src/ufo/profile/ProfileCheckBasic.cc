/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBasic.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBasic> makerProfileCheckBasic_("Basic");

  ProfileCheckBasic::ProfileCheckBasic(const ProfileConsistencyCheckParameters &options,
                                       const ProfileIndices &profileIndices,
                                       ProfileDataHandler &profileDataHandler,
                                       ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckBasic::runCheck()
  {
    oops::Log::debug() << " Basic checks" << std::endl;

    // Set basic check result to true
    result_ = true;

    // Skip this routine if specifically requested
    if (options_.BChecks_Skip.value())
      {
        oops::Log::warning() << "Skipping basic checks" << std::endl;
        return;
      }

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures =
      profileDataHandler_.get<float>(ufo::VariableNames::name_air_pressure);
    // All QC flags are retrieved for the basic checks.
    // (Some might be empty; that is checked before they are used.)
    std::vector <int> &tFlags = profileDataHandler_.get<int>(ufo::VariableNames::name_qc_tFlags);
    std::vector <int> &zFlags = profileDataHandler_.get<int>(ufo::VariableNames::name_qc_zFlags);
    std::vector <int> &uFlags = profileDataHandler_.get<int>(ufo::VariableNames::name_qc_uFlags);

    // Warn and exit if pressures vector is empty
    if (pressures.empty()) {
      result_ = false;
      oops::Log::warning() << "Pressures vector is empty" << std::endl;
      return;
     }

    // Is the number of levels to check OK?
    bool numLevelsToCheckOK = (numLevelsToCheck > 0);

    // Are any levels in the wrong order?
    bool pressOrderOK = true;
    for (int jlev = 0; jlev < numLevelsToCheck - 1; ++jlev) {
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

    oops::Log::debug() << " -> numLevelsToCheckOK: " << numLevelsToCheckOK << std::endl;
    oops::Log::debug() << " -> pressOrderOK: " << pressOrderOK << std::endl;
    oops::Log::debug() << " -> maxPressOK: " << maxPressOK << std::endl;
    oops::Log::debug() << " -> minPressOK: " << minPressOK << std::endl;

    result_ = numLevelsToCheckOK && pressOrderOK && maxPressOK && minPressOK;
    oops::Log::debug() << " -> basicResult: " << result_ << std::endl;

    // If the basic checks are failed, set reject flags
    // This is not done in the OPS sonde consistency checks, but is done in Ops_SondeAverage.inc
    if (options_.flagBasicChecksFail.value() && !result_) {
      for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
        if (!tFlags.empty()) tFlags[jlev] |= ufo::FlagsElem::FinalRejectFlag;
        if (!zFlags.empty()) zFlags[jlev] |= ufo::FlagsElem::FinalRejectFlag;
        if (!uFlags.empty()) uFlags[jlev] |= ufo::FlagsElem::FinalRejectFlag;
      }
    }
  }
}  // namespace ufo

