/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckRH.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckRH> makerProfileCheckRH_("RH");

  ProfileCheckRH::ProfileCheckRH
  (const ProfileConsistencyCheckParameters &options,
   const ProfileIndices &profileIndices,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckRH::runCheck()
  {
    oops::Log::debug() << " Relative humidity check" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::name_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler_.get<float>(ufo::VariableNames::name_obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::name_hofx_air_temperature);
    const std::vector <float> &RHObs =
       profileDataHandler_.get<float>(ufo::VariableNames::name_obs_relative_humidity);
    const std::vector <float> &RHBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::name_hofx_relative_humidity);
    const std::vector <float> &tdObs =
       profileDataHandler_.get<float>(ufo::VariableNames::name_obs_dew_point_temperature);
    const std::vector <int> &tFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_qc_tFlags);
    std::vector <int> &RHFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_qc_RHFlags);
    const std::vector <float> &tObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::name_tObsCorrection);

    std::vector <int> &TotCProfs =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_TotCProfs);
    std::vector <int> &TotHProfs =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_TotHProfs);
    std::vector <int> &TotCFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_TotCFlags);
    std::vector <int> &TotHFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_TotHFlags);
    std::vector <int> &TotLFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_TotLFlags);

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
    float PTrop = 0.0;  // Tropopause pressure level
    int NLowP = 0;  // Number of RH reports above 100 hPa
    float RHDLowP = 0.0;  // Mean RH O-B above 100 hPa
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
          NLowP++;
          RHDLowP = RHDLowP + rh_[NumLev] - rhbk_[NumLev];
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
        float MinRHabove = rh_[jlev];  // Min RH in levels close in pressure to test level
        for (int klev = jlev + 1; klev < NumLev; ++klev) {
          if (Press_[jlev] - Press_[klev] > PressDiffAdjThresh) break;
          MinRHabove = std::min(MinRHabove, rh_[klev]);
        }
        if (MinRHabove < options_.RHCheck_MinRHThresh.value()) {
          FlagH_[jlev] = 2;
          TotCProfs[0]++;
          oops::Log::debug() << " -> Error at top of cloud layer for level " << jlev << std::endl;
          for (int klev = jlev + 1; klev < NumLev; ++klev) {
            if (Press_[jlev] - Press_[klev] > PressDiffAdjThresh) break;
            if (td_[klev] <= td_[jlev - 1]) break;
            FlagH_[klev] = 2;
          }
        }
      }
    }

    if (NLowP > 0) RHDLowP = RHDLowP / static_cast <float> (NLowP);

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
            (Press_[ilev] <= options_.RHCheck_PressInitThresh.value() && RHDLowP > SondeRHHiTol)) {
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

    if (NumHFlags > 0) TotHProfs[0]++;

    for (int n = 0; n < NumCFlags; ++n)
      TotCFlags[0]++;
    for (int n = 0; n < NumHFlags; ++n)
      TotHFlags[0]++;
    for (int n = 0; n < NumLFlags; ++n)
      TotLFlags[0]++;

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
    profileDataHandler_.set(ufo::VariableNames::name_Press, std::move(Press_));
    profileDataHandler_.set(ufo::VariableNames::name_Temp, std::move(Temp_));
    profileDataHandler_.set(ufo::VariableNames::name_rh, std::move(rh_));
    profileDataHandler_.set(ufo::VariableNames::name_td, std::move(td_));
    profileDataHandler_.set(ufo::VariableNames::name_tbk, std::move(tbk_));
    profileDataHandler_.set(ufo::VariableNames::name_rhbk, std::move(rhbk_));
    profileDataHandler_.set(ufo::VariableNames::name_FlagH, std::move(FlagH_));
    profileDataHandler_.set(ufo::VariableNames::name_Indx, std::move(Indx_));
  }
}  // namespace ufo
