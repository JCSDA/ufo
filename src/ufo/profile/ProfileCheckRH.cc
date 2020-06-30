/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckRH.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckRH> makerProfileCheckRH_("RH");

  ProfileCheckRH::ProfileCheckRH
  (const ProfileConsistencyCheckParameters &options,
   const ProfileIndices &profileIndices,
   const ProfileData &profileData,
   ProfileFlags &profileFlags,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator)
  {}

  void ProfileCheckRH::runCheck()
  {
    oops::Log::debug() << " Relative humidity check" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    const std::vector <float> &tObs = profileData_.gettObs();
    const std::vector <float> &tBkg = profileData_.gettBkg();
    const std::vector <float> &RHObs = profileData_.getRHObs();
    const std::vector <float> &RHBkg = profileData_.getRHBkg();
    const std::vector <float> &tdObs = profileData_.gettdObs();
    const std::vector <int> &tFlags = profileFlags_.gettFlags();
    std::vector <int> &RHFlags = profileFlags_.getRHFlags();
    const std::vector <float> &tObsCorrection = profileFlags_.gettObsCorrection();

    if (oops::anyVectorEmpty(pressures, tObs, tBkg, RHObs, RHBkg,
                             tdObs, tFlags, RHFlags, tObsCorrection)) {
      oops::Log::warning() << "At least one vector is empty. "
                           << "Check will not be performed." << std::endl;
      return;
    }
    if (!oops::allVectorsSameSize(pressures, tObs, tBkg, RHObs, RHBkg,
                                  tdObs, tFlags, RHFlags, tObsCorrection)) {
      oops::Log::warning() << "Not all vectors have the same size. "
                           << "Check will not be performed." << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    // Set up arrays and counters

    int NumLev = 0;
    float PTrop = 0.0;
    int Nover100 = 0;
    float RHDover100 = 0.0;
    Press_.assign(numLevelsToCheck, 0.0);
    Temp_.assign(numLevelsToCheck, 0.0);
    rh_.assign(numLevelsToCheck, 0.0);
    td_.assign(numLevelsToCheck, 0.0);
    tbk_.assign(numLevelsToCheck, 0.0);
    rhbk_.assign(numLevelsToCheck, 0.0);
    FlagH_.assign(numLevelsToCheck, 0);
    Indx_.assign(numLevelsToCheck, -1);

    float Tmin = options_.RHCheck_TminInit.value();
    for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
      if (tFlags[jlev] & ufo::FlagsProfile::TropopauseFlag) {
        PTrop = pressures[jlev] * 0.01;
      }
      if (pressures[jlev] > 0.0 &&
          RHObs[jlev] != missingValueFloat) {
        Tmin = std::min(Tmin, tObsFinal[jlev]);
        Indx_[NumLev] = jlev;
        Press_[NumLev] = pressures[jlev] * 0.01;
        Temp_[NumLev] = tObsFinal[jlev];
        rh_[NumLev] = RHObs[jlev];
        td_[NumLev] = tdObs[jlev];
        tbk_[NumLev] = tBkg[jlev];
        rhbk_[NumLev] = RHBkg[jlev];
        if (Press_[NumLev] <= options_.RHCheck_PressInitThresh.value()) {
          Nover100++;
          RHDover100 = RHDover100 + rh_[NumLev] - rhbk_[NumLev];
        }
        NumLev++;
      }
    }
    if (NumLev <= 1) return;

    // Look for errors at the top of cloud layers

    const float RHThresh = options_.RHCheck_RHThresh.value();
    const float PressDiffAdjThresh = options_.RHCheck_PressDiffAdjThresh.value();
    for (int jlev = 1; jlev < NumLev; ++jlev) {
      if (Press_[jlev] < options_.RHCheck_PressThresh.value()) break;
      if ((Press_[0] - Press_[jlev]) < options_.RHCheck_PressDiff0Thresh.value()) continue;
      if (FlagH_[jlev] == 1) continue;
      if ((td_[jlev] - td_[jlev - 1] > options_.RHCheck_tdDiffThresh.value() &&
           rh_[jlev - 1] >= RHThresh) ||
          (td_[jlev] > td_[jlev - 1] &&
           std::min(rh_[jlev - 1], rh_[jlev]) >= RHThresh)) {
        float MinRHabove = rh_[jlev];
        for (int klev = jlev + 1; klev < NumLev; ++klev) {
          if (Press_[jlev] - Press_[klev] > PressDiffAdjThresh) break;
          MinRHabove = std::min(MinRHabove, rh_[klev]);
        }
        if (MinRHabove < options_.RHCheck_MinRHThresh.value()) {
          FlagH_[jlev] = 2;
          profileFlags_.incrementCounterCumul("TotCProfs");
          oops::Log::debug() << " -> Error at top of cloud layer for level " << jlev << std::endl;
          for (int klev = jlev + 1; klev < NumLev; ++klev) {
            if (Press_[jlev] - Press_[klev] > PressDiffAdjThresh) break;
            if (td_[klev] <= td_[jlev - 1]) break;
            FlagH_[klev] = 2;
          }
        }
      }
    }

    if (Nover100 > 0) RHDover100 = RHDover100 / static_cast <float> (Nover100);

    // Simple check for sonde ascent too moist at high levels
    // Start at top and work down

    const float SondeRHHiTol = options_.RHCheck_SondeRHHiTol.value();

    int NumLFlags = 0;
    int NumHFlags = 0;
    int NumCFlags = 0;
    if ((PTrop != 0.0 && Press_[NumLev - 1] <= PTrop) ||
        Tmin < options_.RHCheck_TminThresh.value()) {
      for (int ilev = NumLev - 1; ilev >= 0; ilev--) {
        if (rh_[ilev] > rhbk_[ilev] + SondeRHHiTol ||
            (Press_[ilev] <= 100.0 && RHDover100 > SondeRHHiTol)) {
          if (FlagH_[ilev] == 0) {
            oops::Log::debug() << " -> Sonde ascent too moist for level " << ilev << std::endl;
            FlagH_[ilev] = 1;
            if (Temp_[ilev] >= options_.RHCheck_TempThresh.value()) NumLFlags = NumLFlags + 1;
          }
        } else {
          break;
        }
      }
    }

    // Set level flags
    for (int i = 0; i < NumLev; ++i) {
      if (FlagH_[i] == 1) ++NumHFlags;
      if (FlagH_[i] == 2) ++NumCFlags;
    }

    if (NumHFlags > 0) profileFlags_.incrementCounterCumul("TotHProfs");

    for (int n = 0; n < NumCFlags; ++n)
      profileFlags_.incrementCounterCumul("TotCFlags");
    for (int n = 0; n < NumHFlags; ++n)
      profileFlags_.incrementCounterCumul("TotHFlags");
    for (int n = 0; n < NumLFlags; ++n)
      profileFlags_.incrementCounterCumul("TotLFlags");

    if (NumCFlags + NumHFlags > 0) {
      for (int jlev = 0; jlev < NumLev; ++jlev) {
        if (FlagH_[jlev] > 0) {
          int ilev = Indx_[jlev];
          RHFlags[ilev] |= ufo::FlagsProfile::InterpolationFlag;
          RHFlags[ilev] |= ufo::FlagsElem::FinalRejectFlag;
        }
      }
    }
  }

  void ProfileCheckRH::fillValidator()
  {
    profileCheckValidator_.setRHFlags(profileFlags_.getuFlags());
    profileCheckValidator_.setTotCProfs(profileFlags_.getCounter("TotCProfs"));
    profileCheckValidator_.setTotHProfs(profileFlags_.getCounter("TotHProfs"));
    profileCheckValidator_.setTotCFlags(profileFlags_.getCounter("TotCFlags"));
    profileCheckValidator_.setTotHFlags(profileFlags_.getCounter("TotHFlags"));
    profileCheckValidator_.setTotLFlags(profileFlags_.getCounter("TotLFlags"));
    profileCheckValidator_.setPress(Press_);
    profileCheckValidator_.setTemp(Temp_);
    profileCheckValidator_.setrh(rh_);
    profileCheckValidator_.settd(td_);
    profileCheckValidator_.settbk(tbk_);
    profileCheckValidator_.setrhbk(rhbk_);
    profileCheckValidator_.setFlagH(FlagH_);
    profileCheckValidator_.setIndx(Indx_);
  }
}  // namespace ufo
