/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckUInterp.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckUInterp> makerProfileCheckUInterp_("UInterp");

  ProfileCheckUInterp::ProfileCheckUInterp
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckUInterp::runCheck()
  {
    oops::Log::debug() << " U interpolation check" << std::endl;

    const int numProfileLevels = profileDataHandler_.getNumProfileLevels();
    const std::vector <float> &pressures =
      profileDataHandler_.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &uObs =
      profileDataHandler_.get<float>(ufo::VariableNames::obs_eastward_wind);
    const std::vector <float> &vObs =
      profileDataHandler_.get<float>(ufo::VariableNames::obs_northward_wind);
    std::vector <int> &uFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &NumSamePErrObs =
      profileDataHandler_.get<int>(ufo::VariableNames::counter_NumSamePErrObs);
    std::vector <int> &NumInterpErrObs =
      profileDataHandler_.get<int>(ufo::VariableNames::counter_NumInterpErrObs);

    if (!oops::allVectorsSameNonZeroSize(pressures, uObs, vObs, uFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(pressures, uObs, vObs, uFlags)
                           << std::endl;
      return;
    }

    calcStdLevelsUV(numProfileLevels, pressures, uObs, vObs, uFlags);

    LevErrors_.assign(numProfileLevels, -1);
    uInterp_.assign(numProfileLevels, 0.0);
    vInterp_.assign(numProfileLevels, 0.0);

    int NumErrors = 0;

    // Check levels with identical pressures
    int jlevprev = -1;
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (uObs[jlev] != missingValueFloat && vObs[jlev] != missingValueFloat) {
        if (jlevprev != -1) {
          if (pressures[jlev] == pressures[jlevprev] &&
              (uObs[jlev] != uObs[jlevprev] || vObs[jlev] != vObs[jlevprev])) {
            float VectDiffSq = std::pow(uObs[jlev] - uObs[jlevprev], 2) +
              std::pow(vObs[jlev] - vObs[jlevprev], 2);
            if (VectDiffSq > options_.UICheck_TInterpIdenticalPTolSq.value()) {
              NumErrors++;
              uFlags[jlevprev] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
              uFlags[jlev]     |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
              oops::Log::debug() << " -> Wind speed interpolation check: identical P for "
                                 << "levels " << jlevprev << " and " << jlev << std::endl;
              oops::Log::debug() << " -> Level " << jlevprev << ": "
                                 << "P = " << pressures[jlevprev] * 0.01 << "hPa, uObs = "
                                 << uObs[jlevprev] << "ms^-1, vObs = "
                                 << vObs[jlevprev] << "ms^-1" << std::endl;
              oops::Log::debug() << " -> Level " << jlev << ": "
                                 << "P = " << pressures[jlev] * 0.01 << "hPa, uObs = "
                                 << uObs[jlev] << "ms^-1, vObs = "
                                 << vObs[jlev] << "ms^-1" << std::endl;
              oops::Log::debug() << " -> VectDiffSq = " << VectDiffSq << "m^2s^-2" << std::endl;
            }
          }
        }
        jlevprev = jlev;
      }
    }

    if (NumErrors > 0) NumSamePErrObs[0]++;

    if (NumSig_ < std::max(3, NumStd_ / 2)) return;  // Too few sig levels for reliable check

    // Interpolation check
    NumErrors = 0;
    for (int jlevStd = 0; jlevStd < NumStd_; ++jlevStd) {
      if (SigBelow_[jlevStd] == -1 || SigAbove_[jlevStd] == -1) continue;
      int jlev = StdLev_[jlevStd];  // Standard level
      int SigB = SigBelow_[jlevStd];
      int SigA = SigAbove_[jlevStd];
      float PStd = pressures[jlev];
      // BigGap - see 6.3.2.2.2 of the Guide on the Global Data-Processing System
      float BigGap = options_.UICheck_BigGapLowP.value();
      const std::vector <float> BigGaps = options_.UICheck_BigGaps.value();
      const std::vector <float> BigGapsPThresh = options_.UICheck_BigGapsPThresh.value();
      for (size_t bgidx = 0; bgidx < BigGapsPThresh.size(); ++bgidx) {
        if (PStd > BigGapsPThresh[bgidx]) {
          BigGap = BigGaps[bgidx];
          break;
        }
      }

      if (pressures[SigB] - PStd > BigGap ||
          PStd - pressures[SigA] > BigGap) {
        uInterp_[jlev] = missingValueFloat;
        vInterp_[jlev] = missingValueFloat;
        continue;
      }

      if (LogP_[SigB] == LogP_[SigA]) continue;

      float Ratio = (LogP_[jlev] - LogP_[SigB]) / (LogP_[SigA] - LogP_[SigB]);  // eqn 3.3a
      uInterp_[jlev] = uObs[SigB] + (uObs[SigA] - uObs[SigB]) * Ratio;  // eqn 3.3b
      vInterp_[jlev] = vObs[SigB] + (vObs[SigA] - vObs[SigB]) * Ratio;  // eqn 3.3b

      // Vector wind difference > UInterpTol m/s?
      float VectDiffSq = std::pow(uObs[jlev] - uInterp_[jlev], 2) +
        std::pow(vObs[jlev] - vInterp_[jlev], 2);
      if (VectDiffSq > options_.UICheck_TInterpTolSq.value()) {
        NumErrors++;
        LevErrors_[jlev]++;
        LevErrors_[SigB]++;
        LevErrors_[SigA]++;
        // Simplest form of flagging
        uFlags[jlev] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
        uFlags[SigB] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
        uFlags[SigA] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;

        oops::Log::debug() << " -> Failed wind speed interpolation check for levels " << jlev
                           << " (central), " << SigB << " (lower) and "
                           << SigA << " (upper)" << std::endl;
        oops::Log::debug() << " -> Level " << jlev << ": "
                           << "P = " << pressures[jlev] * 0.01 << "hPa, uObs = "
                           << uObs[jlev] << "ms^-1, vObs = "
                           << vObs[jlev] << "ms^-1, uInterp = " << uInterp_[jlev]
                           << "ms^-1, vInterp = " << vInterp_[jlev] << "ms^-1" << std::endl;
        oops::Log::debug() << " -> VectDiffSq = " << VectDiffSq << "m^2s^-2" << std::endl;
      }
    }

    if (NumErrors > 0) NumInterpErrObs[0]++;
  }

  void ProfileCheckUInterp::fillValidator()
  {
    profileDataHandler_.set(ufo::VariableNames::StdLev, std::move(StdLev_));
    profileDataHandler_.set(ufo::VariableNames::SigAbove, std::move(SigAbove_));
    profileDataHandler_.set(ufo::VariableNames::SigBelow, std::move(SigBelow_));
    profileDataHandler_.set(ufo::VariableNames::LevErrors, std::move(LevErrors_));
    profileDataHandler_.set(ufo::VariableNames::uInterp, std::move(uInterp_));
    profileDataHandler_.set(ufo::VariableNames::vInterp, std::move(vInterp_));
    profileDataHandler_.set(ufo::VariableNames::LogP, std::move(LogP_));
    std::vector <int> NumStd(profileDataHandler_.getNumProfileLevels(), std::move(NumStd_));
    std::vector <int> NumSig(profileDataHandler_.getNumProfileLevels(), std::move(NumSig_));
    profileDataHandler_.set(ufo::VariableNames::NumStd, std::move(NumStd));
    profileDataHandler_.set(ufo::VariableNames::NumSig, std::move(NumSig));
  }
}  // namespace ufo
