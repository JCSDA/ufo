/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBasic.h"

namespace ufo {
  ProfileCheckBasic::ProfileCheckBasic(const ProfileConsistencyCheckParameters &options,
                                       const ProfileIndices &profileIndices,
                                       const ProfileData &profileData,
                                       ProfileFlags &profileFlags,
                                       ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator)
  {}

  void ProfileCheckBasic::runCheck()
  {
    oops::Log::debug() << " Basic checks" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    std::vector <int> &tFlags = profileFlags_.gettFlags();

    // Is the number of levels to check OK?
    // (Original OPS condition is (numLevelsToCheck != 0 && numLevelsToCheck != IMDI_))
    bool numLevelsToCheckOK = (numLevelsToCheck > 0);

    // Are any levels in the wrong order?
    bool pressOrderOK = true;
    for (int jlev = 0; jlev < numLevelsToCheck - 1; ++jlev) {
      pressOrderOK = pressOrderOK && (pressures[jlev] >= pressures[jlev + 1]);
      if (!pressOrderOK) break;
    }

    // Is the first level > maximum value pressure?
    bool maxPressOK = (pressures.size() > 0 ?
                       pressures.front() <= options_.BChecks_maxValidP.value() :
                       false);

    // Is the last level < minimum value pressure?
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
        tFlags[jlev] |= ufo::FlagsElem::FinalRejectFlag;
      }
    }
  }

  void ProfileCheckBasic::fillValidator()
  {
    profileCheckValidator_.settFlags(profileFlags_.gettFlags());
    profileCheckValidator_.setNumAnyErrors(profileFlags_.getCounter("NumAnyErrors"));
  }
}  // namespace ufo

