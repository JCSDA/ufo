/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckSign.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSign> makerProfileCheckSign_("Sign");

  ProfileCheckSign::ProfileCheckSign(const ProfileConsistencyCheckParameters &options,
                                     const ProfileIndices &profileIndices,
                                     ProfileDataHandler &profileDataHandler,
                                     ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckSign::runCheck()
  {
    oops::Log::debug() << " Sign check/correction" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();

    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::name_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler_.get<float>(ufo::VariableNames::name_obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::name_hofx_air_temperature);
    const std::vector <float> &PstarBackgr =
       profileDataHandler_.get<float>(ufo::VariableNames::name_PstarBackgr);
    std::vector <int> &tFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_qc_tFlags);
    std::vector <int> &NumAnyErrors =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_NumAnyErrors);
    std::vector <int> &NumSignChange =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_NumSignChange);
    std::vector <float> &tObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::name_tObsCorrection);
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
          NumAnyErrors[0]++;
          NumSignChange[0]++;

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
}  // namespace ufo

