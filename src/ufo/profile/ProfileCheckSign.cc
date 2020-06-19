/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckSign.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSign> makerProfileCheckSign_("Sign");

  ProfileCheckSign::ProfileCheckSign(const ProfileConsistencyCheckParameters &options,
                                     const ProfileIndices &profileIndices,
                                     const ProfileData &profileData,
                                     ProfileFlags &profileFlags,
                                     ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator)
  {}

  void ProfileCheckSign::runCheck()
  {
    oops::Log::debug() << " Sign check/correction" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    const std::vector <float> &tObs = profileData_.gettObs();
    const std::vector <float> &tBkg = profileData_.gettBkg();
    const std::vector <float> &PstarBackgr = profileData_.getPstarBackgr();
    std::vector <int> &tFlags = profileFlags_.gettFlags();
    std::vector <float> &tObsCorrection =
      profileFlags_.gettObsCorrection();  // Potentially modified here

    if (oops::anyVectorEmpty(pressures, tObs, tBkg, PstarBackgr, tFlags, tObsCorrection)) {
      oops::Log::warning() << "At least one vector is empty. "
                           << "Check will not be performed." << std::endl;
      return;
    }
    if (!oops::allVectorsSameSize(pressures, tObs, tBkg, PstarBackgr, tFlags, tObsCorrection)) {
      oops::Log::warning() << "Not all vectors have the same size. "
                           << "Check will not be performed." << std::endl;
      return;
    }

    for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
      if (tFlags[jlev] & ufo::FlagsElem::FinalRejectFlag) continue;  // Ignore this level
      if (pressures[jlev] <= PstarBackgr[jlev] - options_.SCheck_PstarThresh.value() &&
          tObs[jlev] != missingValueFloat &&
          std::abs(tObs[jlev] - tBkg[jlev]) >= options_.SCheck_tObstBkgThresh.value()) {
        // Change sign of tObs in C and compare to tBkg (also in C)
        if (std::abs(2.0 * ufo::Constants::t0c - tObs[jlev] - tBkg[jlev]) <
            options_.SCheck_ProfileSignTol.value()) {
          profileFlags_.incrementCounter("NumAnyErrors");
          profileFlags_.incrementCounterCumul("NumSignChange");

          tFlags[jlev] |= ufo::FlagsElem::DataCorrectFlag;

          oops::Log::debug() << " -> Failed sign check for level " << jlev << std::endl;
          oops::Log::debug() << " -> P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << ufo::Constants::t0c - tObs[jlev] << "C, tBkg = "
                             << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;

          if (options_.SCheck_CorrectT.value()) {
            // Corrected T is 2 * t0c - T (all quantities in K).
            // The correction is 2 * (t0c - T).
            tObsCorrection[jlev] = 2.0 * (ufo::Constants::t0c - tObs[jlev]);

            oops::Log::debug() << " -> Uncorrected tObs: " << tObs[jlev] << "C" << std::endl;
            oops::Log::debug() << "    tObs correction: "
                               << tObsCorrection[jlev] << "C" << std::endl;
            oops::Log::debug() << "    Corrected tObs: "
                               << tObs[jlev] + tObsCorrection[jlev] << "C" << std::endl;
          } else {
            // Observation is rejected
            tFlags[jlev] |= ufo::FlagsElem::FinalRejectFlag;
          }
        } else if (pressures[jlev] > options_.SCheck_PrintLargeTThresh.value()) {
          // Print out information on other large T differences
          oops::Log::debug() << " -> Passed test but have large T difference for level "
                             << jlev << ": "
                             << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << tObs[jlev] - ufo::Constants::t0c << "C, tBkg = "
                             << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
        }
      }
    }
  }

  void ProfileCheckSign::fillValidator()
  {
    profileCheckValidator_.settFlags(profileFlags_.gettFlags());
    profileCheckValidator_.setNumAnyErrors(profileFlags_.getCounter("NumAnyErrors"));
    profileCheckValidator_.setNumSignChange(profileFlags_.getCounter("NumSignChange"));
  }
}  // namespace ufo

