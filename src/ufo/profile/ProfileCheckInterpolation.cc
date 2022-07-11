/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckInterpolation.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckInterpolation>
  makerProfileCheckInterpolation_("Interpolation");

  ProfileCheckInterpolation::ProfileCheckInterpolation
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckInterpolation::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Interpolation check" << std::endl;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler.get<float>(ufo::VariableNames::hofx_air_temperature);
    std::vector <int> &tFlags =
       profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &NumAnyErrors =
       profileDataHandler.get<int>(ufo::VariableNames::counter_NumAnyErrors);
    std::vector <int> &NumInterpErrors =
       profileDataHandler.get<int>(ufo::VariableNames::counter_NumInterpErrors);
    std::vector <int> &NumInterpErrObs =
       profileDataHandler.get<int>(ufo::VariableNames::counter_NumInterpErrObs);
    const std::vector <float> &tObsCorrection =
       profileDataHandler.get<float>(ufo::VariableNames::obscorrection_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg, tFlags,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg, tFlags,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    calcStdLevels(numProfileLevels, pressures, tObsFinal, tFlags);

    LevErrors_.assign(numProfileLevels, -1);
    tInterp_.assign(numProfileLevels, missingValueFloat);

    int NumErrors = 0;

    for (int jlevstd = 0; jlevstd < NumStd_; ++jlevstd) {
      int jlev = StdLev_[jlevstd];  // Standard level

      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag) continue;
      int SigB = SigBelow_[jlevstd];
      int SigA = SigAbove_[jlevstd];
      float PStd = pressures[jlev];
      int IPStd = std::round(pressures[jlev] * 0.01);  // Pressure rounded to nearest hPa

      /// BigGap - see 6.3.2.2.2 of the Guide on the Global Data-Processing System.
      /// Reduced to 50 hPa for standard levels at 150 and 100 hPa
      float BigGap = options_.ICheck_BigGapInit.value();
      for (int i = 0; i < StandardLevels_.size(); ++i) {
        if (StandardLevels_[i] <= IPStd) {
          BigGap = BigGaps_[i] * 100.0;  // hPa -> Pa
          break;
        }
      }

      if (NumSig_ < std::max(3, NumStd_ / 2)) continue;  // Too few sig levs for reliable check

      if (SigB == -1 || SigA == -1) continue;

      if (pressures[SigB] - PStd > BigGap ||
          PStd - pressures[SigA] > BigGap ||
          LogP_[SigB] == LogP_[SigA]) continue;

      float Ratio = (LogP_[jlev] - LogP_[SigB]) /
        (LogP_[SigA] - LogP_[SigB]);  // eqn 3.3a

      tInterp_[jlev] = tObsFinal[SigB] + (tObsFinal[SigA] - tObsFinal[SigB]) * Ratio;  // eqn 3.3b

      // Temperature difference > TInterpTol*TolRelax degrees?
      float TolRelax = 1.0;
      if (PStd < options_.ICheck_TolRelaxPThresh.value())
        TolRelax = options_.ICheck_TolRelax.value();
      if (std::abs(tObsFinal[jlev] - tInterp_[jlev]) >
          options_.ICheck_TInterpTol.value() * TolRelax) {
        NumAnyErrors[0]++;
        NumInterpErrors[0]++;
        NumErrors++;

        // Simplest form of flagging - sig or std flags may be unset in other routines
        tFlags[jlev] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
        tFlags[SigB] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
        tFlags[SigA] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;

        LevErrors_[jlev]++;
        LevErrors_[SigB]++;
        LevErrors_[SigA]++;

        oops::Log::debug() << " -> Failed interpolation check for levels " << jlev
                           << " (central), " << SigB << " (lower) and "
                           << SigA << " (upper)" << std::endl;
        oops::Log::debug() << " -> Level " << jlev << ": "
                           << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                           << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c << "C, "
                           << "tInterp = " << tInterp_[jlev] - ufo::Constants::t0c
                           << "C, tInterp - tObs = " << tInterp_[jlev] - tObsFinal[jlev]
                           << std::endl;
        oops::Log::debug() << " -> Level " << SigB << ": "
                           << "P = " << pressures[SigB] * 0.01 << "hPa, tObs = "
                           << tObsFinal[SigB] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[SigB] - ufo::Constants::t0c << "C" << std::endl;
        oops::Log::debug() << " -> Level " << SigA << ": "
                           << "P = " << pressures[SigA] * 0.01 << "hPa, tObs = "
                           << tObsFinal[SigA] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[SigA] - ufo::Constants::t0c << "C" << std::endl;
      }
    }
    if (NumErrors > 0) NumInterpErrObs[0]++;
  }

  void ProfileCheckInterpolation::fillValidationData(ProfileDataHandler &profileDataHandler)
  {
    profileDataHandler.set(ufo::VariableNames::StdLev, std::move(StdLev_));
    profileDataHandler.set(ufo::VariableNames::SigAbove, std::move(SigAbove_));
    profileDataHandler.set(ufo::VariableNames::SigBelow, std::move(SigBelow_));
    profileDataHandler.set(ufo::VariableNames::IndStd, std::move(IndStd_));
    profileDataHandler.set(ufo::VariableNames::LevErrors, std::move(LevErrors_));
    profileDataHandler.set(ufo::VariableNames::tInterp, std::move(tInterp_));
    profileDataHandler.set(ufo::VariableNames::LogP, std::move(LogP_));
    std::vector <int> NumStd(profileDataHandler.getNumProfileLevels(), std::move(NumStd_));
    std::vector <int> NumSig(profileDataHandler.getNumProfileLevels(), std::move(NumSig_));
    profileDataHandler.set(ufo::VariableNames::NumStd, std::move(NumStd));
    profileDataHandler.set(ufo::VariableNames::NumSig, std::move(NumSig));
  }
}  // namespace ufo

