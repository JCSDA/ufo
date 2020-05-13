/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckUnstableLayer.h"

namespace ufo {
  ProfileCheckUnstableLayer::ProfileCheckUnstableLayer
  (const ProfileConsistencyCheckParameters &options,
   const ProfileIndices &profileIndices,
   const ProfileData &profileData,
   ProfileFlags &profileFlags,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator)
  {}

  void ProfileCheckUnstableLayer::runCheck()
  {
    oops::Log::debug() << " Unstable layer/superadiabat check" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    const std::vector <float> &tObs = profileData_.gettObs();
    const std::vector <float> &tBkg = profileData_.gettBkg();
    std::vector <int> &tFlags = profileFlags_.gettFlags();
    const std::vector <float> &tObsCorrection = profileFlags_.gettObsCorrection();
    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    PBottom_ = 0.0;

    int jlevprev = 0;
    for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
      if (tFlags[jlev] & ufo::FlagsElem::FinalRejectFlag) continue;  // Ignore this level
      if (tObsFinal[jlev] != missingValueFloat &&
          pressures[jlev] > options_.ULCheck_MinP.value()) {
        if (PBottom_ == 0.0) {
          PBottom_ = pressures[jlev];
        } else {
          //  Temperature calculated adiabatically from prev level
          const float Tadiabat = tObsFinal[jlevprev] *
            std::pow(pressures[jlev] / pressures[jlevprev], ufo::Constants::rd_over_cp);
          if (tObsFinal[jlev] - Tadiabat <= options_.ULCheck_SuperadiabatTol.value() &&
              pressures[jlevprev] <= PBottom_ - options_.ULCheck_PBThresh.value()) {
            profileFlags_.incrementCounter("NumAnyErrors");
            profileFlags_.incrementCounterCumul("NumSuperadiabat");
            tFlags[jlevprev] |= ufo::FlagsProfile::SuperadiabatFlag;
            tFlags[jlev]     |= ufo::FlagsProfile::SuperadiabatFlag;

            oops::Log::debug() << " -> Failed unstable layer/superadiabat check for levels "
                               << jlevprev << " and " << jlev << std::endl;
            oops::Log::debug() << " -> Tadiabat = " << Tadiabat - ufo::Constants::t0c << "C"
                               << std::endl;
            oops::Log::debug() << " -> Level " << jlevprev << ": "
                               << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlevprev] - ufo::Constants::t0c << "C, "
                               << "tObs - Tadiabat = " << tObsFinal[jlev] - Tadiabat
                               << "C" << std::endl;
            oops::Log::debug() << " -> Level " << jlev << ": "
                               << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
          }
        }
        jlevprev = jlev;
      }
    }
  }

  void ProfileCheckUnstableLayer::fillValidator()
  {
    profileCheckValidator_.settFlags(profileFlags_.gettFlags());
    profileCheckValidator_.setNumAnyErrors(profileFlags_.getCounter("NumAnyErrors"));
    profileCheckValidator_.setNumSuperadiabat(profileFlags_.getCounter("NumSuperadiabat"));
    profileCheckValidator_.setPBottom(PBottom_);
  }
}  // namespace ufo
