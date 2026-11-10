/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckInterpolation.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckInterpolation>
  makerProfileCheckInterpolation_("Interpolation");

  ProfileCheckInterpolation::ProfileCheckInterpolation
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckInterpolation::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Interpolation check" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_surface_level_t,
       ufo::ProfileVariableNames::diagflags_standard_level_t,
       ufo::ProfileVariableNames::diagflags_interpolation_t,
       ufo::ProfileVariableNames::diagflags_final_reject_t};
    const std::vector<std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumAnyErrors,
       ufo::ProfileVariableNames::counter_NumInterpErrors,
       ufo::ProfileVariableNames::counter_NumInterpErrObs};
    const std::vector<std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::obscorrection_air_temperature};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
    for (const auto &variableInt : variableNamesInt) {
      profileDataHandler.get<int>(variableInt);
    }
    for (const auto &variableFloat : variableNamesFloat) {
      profileDataHandler.get<float>(variableFloat);
    }
  }

  void ProfileCheckInterpolation::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsTSurfaceLevel =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_t);
    std::vector <bool> &diagFlagsTStandardLevel =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_t);
    std::vector <bool> &diagFlagsTInterp =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_t);
    std::vector <bool> &diagFlagsTFinalReject =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <int> &NumAnyErrors =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumAnyErrors);
    std::vector <int> &NumInterpErrors =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumInterpErrors);
    std::vector <int> &NumInterpErrObs =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumInterpErrObs);
    const std::vector <float> &tObsCorrection =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obscorrection_air_temperature);

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg,
                                         diagFlagsTSurfaceLevel,
                                         diagFlagsTStandardLevel,
                                         diagFlagsTInterp,
                                         diagFlagsTFinalReject,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg,
                                                    diagFlagsTSurfaceLevel,
                                                    diagFlagsTStandardLevel,
                                                    diagFlagsTInterp,
                                                    diagFlagsTFinalReject,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    calcStdLevels(numProfileLevels, pressures, tObsFinal,
                  diagFlagsTFinalReject,
                  diagFlagsTSurfaceLevel,
                  diagFlagsTStandardLevel);

    LevErrors_.assign(numProfileLevels, -1);
    tInterp_.assign(numProfileLevels, missingValueFloat);

    int NumErrors = 0;

    for (int jlevstd = 0; jlevstd < NumStd_; ++jlevstd) {
      int jlev = StdLev_[jlevstd];  // Standard level

      if (diagFlagsTSurfaceLevel[jlev]) continue;
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
        diagFlagsTInterp[jlev] = true;
        diagFlagsTInterp[SigB] = true;
        diagFlagsTInterp[SigA] = true;

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
    profileDataHandler.set(ufo::ProfileVariableNames::StdLev, std::move(StdLev_));
    profileDataHandler.set(ufo::ProfileVariableNames::SigAbove, std::move(SigAbove_));
    profileDataHandler.set(ufo::ProfileVariableNames::SigBelow, std::move(SigBelow_));
    profileDataHandler.set(ufo::ProfileVariableNames::IndStd, std::move(IndStd_));
    profileDataHandler.set(ufo::ProfileVariableNames::LevErrors, std::move(LevErrors_));
    profileDataHandler.set(ufo::ProfileVariableNames::tInterp, std::move(tInterp_));
    profileDataHandler.set(ufo::ProfileVariableNames::LogP, std::move(LogP_));
    std::vector <int> NumStd(profileDataHandler.getNumProfileLevels(), std::move(NumStd_));
    std::vector <int> NumSig(profileDataHandler.getNumProfileLevels(), std::move(NumSig_));
    profileDataHandler.set(ufo::ProfileVariableNames::NumStd, std::move(NumStd));
    profileDataHandler.set(ufo::ProfileVariableNames::NumSig, std::move(NumSig));
  }
}  // namespace ufo

