/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckSamePDiffT.h"

namespace ufo {
  ProfileCheckSamePDiffT::ProfileCheckSamePDiffT(const ProfileConsistencyCheckParameters &options,
                                                 const ProfileIndices &profileIndices,
                                                 const ProfileData &profileData,
                                                 ProfileFlags &profileFlags,
                                                 ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator)
  {}

  void ProfileCheckSamePDiffT::runCheck()
  {
    oops::Log::debug() << " Test for same pressure and different temperature" << std::endl;
    int jlevprev = -1;
    int NumErrors = 0;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    const std::vector <float> &tObs = profileData_.gettObs();
    const std::vector <float> &tBkg = profileData_.gettBkg();
    std::vector <int> &tFlags = profileFlags_.gettFlags();
    const std::vector <float> &tObsCorrection = profileFlags_.gettObsCorrection();
    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
      if (tObs[jlev] == missingValueFloat) continue;

      if (jlevprev == -1) {
        jlevprev = jlev;
        continue;
      }

      if (pressures[jlev] == pressures[jlevprev]) {
        int jlevuse = jlevprev;
        if (std::abs(tObsFinal[jlev] - tObsFinal[jlevprev]) > options_.SPDTCheck_TThresh.value()) {
          NumErrors++;
          profileFlags_.incrementCounter("NumAnyErrors");

          // Choose which level to flag
          if (std::abs(tObsFinal[jlev] - tBkg[jlev]) <=
              std::abs(tObsFinal[jlevprev] - tBkg[jlevprev])) {
            tFlags[jlevprev] |= ufo::FlagsElem::FinalRejectFlag;
            tFlags[jlev]     |= ufo::FlagsProfile::InterpolationFlag;
            jlevuse = jlev;
          } else {
            tFlags[jlevprev] |= ufo::FlagsProfile::InterpolationFlag;
            tFlags[jlev]     |= ufo::FlagsElem::FinalRejectFlag;
          }

          oops::Log::debug() << " -> Failed same P/different T check for levels "
                             << jlevprev << " and " << jlev << std::endl;
          oops::Log::debug() << " -> Level " << jlevprev << ": "
                             << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlevprev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> Level " << jlev << ": "
                             << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> tObs difference: " << tObsFinal[jlev] - tObsFinal[jlevprev]
                             << std::endl;
          oops::Log::debug() << " -> Use level " << jlevuse << std::endl;
        }
        jlevprev = jlevuse;
      } else {  // Distinct pressures
        jlevprev = jlev;
      }
    }
    if (NumErrors > 0) profileFlags_.incrementCounterCumul("NumSamePErrObs");
  }

  void ProfileCheckSamePDiffT::fillValidator()
  {
    profileCheckValidator_.settFlags(profileFlags_.gettFlags());
    profileCheckValidator_.setNumAnyErrors(profileFlags_.getCounter("NumAnyErrors"));
    profileCheckValidator_.setNumSamePErrObs(profileFlags_.getCounter("NumSamePErrObs"));
  }
}  // namespace ufo

