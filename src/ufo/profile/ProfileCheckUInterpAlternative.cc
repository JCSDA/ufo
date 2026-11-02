/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckUInterpAlternative.h"
#include "ufo/profile/ProfileVariableNames.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckUInterpAlternative>
  makerProfileCheckUInterpAlternative_("UInterpAlternative");

  ProfileCheckUInterpAlternative::ProfileCheckUInterpAlternative
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckUInterpAlternative::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Alternative U interpolation check" << std::endl;

    // Produce vector of profiles containing data for the alternative U interpolation check.
    std::vector <std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumSamePErrObs,
       ufo::ProfileVariableNames::counter_NumInterpErrObs,
       ufo::ProfileVariableNames::extended_obs_space};
    std::vector <std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::obs_air_pressure,
       ufo::ProfileVariableNames::obs_eastward_wind,
       ufo::ProfileVariableNames::obs_northward_wind};
    std::vector<std::string> variableNamesBool = {
      ufo::ProfileVariableNames::diagflags_interpolation_u,
      ufo::ProfileVariableNames::diagflags_standard_level_u};
    if (options_.compareWithOPS.value()) {
      variableNamesInt.insert(variableNamesInt.end(),
                              {ufo::ProfileVariableNames::StdLev,
                                  ufo::ProfileVariableNames::SigAbove,
                                  ufo::ProfileVariableNames::SigBelow,
                                  ufo::ProfileVariableNames::LevErrors,
                                  ufo::ProfileVariableNames::NumStd,
                                  ufo::ProfileVariableNames::NumSig});
      variableNamesFloat.insert(variableNamesFloat.end(),
                                {ufo::ProfileVariableNames::uInterp,
                                    ufo::ProfileVariableNames::vInterp,
                                    ufo::ProfileVariableNames::LogP});
    }

    // Ensure variables are added to the ObsSpace on ranks with zero obs.
    for (const auto &variable : variableNamesInt) {
      profileDataHandler.get<int>(variable);
    }
    for (const auto &variable : variableNamesFloat) {
      profileDataHandler.get<float>(variable);
    }
    for (const auto &variable : variableNamesBool) {
      profileDataHandler.get<bool>(variable);
    }
    profileDataHandler.get<bool>(
        ufo::ProfileVariableNames::diagflags_final_reject_report);
    profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumAnyErrors);

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
       variableNamesBool,
       {});

    // Run alternative U interpolation check on each original profile.
    const size_t nprofs = profileDataHandler.getObsdb().nrecs();
    for (size_t jprof = 0; jprof < nprofs; ++jprof) {
      oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;
      auto& profile = profiles[jprof];
      // Check whether this profile is in the original ObsSpace or has been averaged
      // onto model levels. If the former, proceed with the check.
      // If the ObsSpace has not been extended then all profiles are by default in
      // the original ObsSpace.
      const auto &extended_obs_space =
        profile.get<int>(ufo::ProfileVariableNames::extended_obs_space);
      if (extended_obs_space.empty() ||
          std::find(extended_obs_space.begin(), extended_obs_space.end(), 0) !=
          extended_obs_space.end())
        runCheckOnProfile(profile);
      // Fill validation information if required.
      if (options_.compareWithOPS.value())
        fillValidationData(profile);
    }

    // Update data handler with profile information.
    oops::Log::debug() << " Updating data handler" << std::endl;
    profileDataHandler.updateAllProfiles(profiles);
  }

  void ProfileCheckUInterpAlternative::runCheckOnProfile(ProfileDataHolder &profile)
  {
    const int numProfileLevels = profile.getNumProfileLevels();
    const std::vector <float> &pressures =
      profile.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <float> &uObs =
      profile.get<float>(ufo::ProfileVariableNames::obs_eastward_wind);
    const std::vector <float> &vObs =
      profile.get<float>(ufo::ProfileVariableNames::obs_northward_wind);
    std::vector <bool> &diagFlagsUInterp =
      profile.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_u);
    const std::vector <bool> &diagFlagsUStdLev =
      profile.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_u);
    std::vector <int> &NumSamePErrObs =
      profile.get<int>(ufo::ProfileVariableNames::counter_NumSamePErrObs);
    std::vector <int> &NumInterpErrObs =
      profile.get<int>(ufo::ProfileVariableNames::counter_NumInterpErrObs);

    if (pressures.empty() ||
        pressures.front() > options_.BChecks_maxValidP.value() ||
        pressures.back() < options_.BChecks_minValidP.value()) {
      return;
    }

    if (!oops::allVectorsSameNonZeroSize(pressures, uObs, vObs,
                                         diagFlagsUInterp,
                                         diagFlagsUStdLev)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, uObs, vObs,
                                                    diagFlagsUInterp,
                                                    diagFlagsUStdLev)
                         << std::endl;
      return;
    }

    calcStdLevelsUV(numProfileLevels, pressures, uObs, vObs,
                    diagFlagsUStdLev);

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
              diagFlagsUInterp[jlevprev] = true;
              diagFlagsUInterp[jlev] = true;
              oops::Log::debug() << " -> Wind speed interpolation check: identical P "
                                 << "and significantly different wind speed magnitude for "
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
        diagFlagsUInterp[jlev] = true;
        diagFlagsUInterp[SigB] = true;
        diagFlagsUInterp[SigA] = true;

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

  void ProfileCheckUInterpAlternative::fillValidationData(ProfileDataHolder &profile)
  {
    profile.set(ufo::ProfileVariableNames::StdLev, std::move(StdLev_));
    profile.set(ufo::ProfileVariableNames::SigAbove, std::move(SigAbove_));
    profile.set(ufo::ProfileVariableNames::SigBelow, std::move(SigBelow_));
    profile.set(ufo::ProfileVariableNames::LevErrors, std::move(LevErrors_));
    profile.set(ufo::ProfileVariableNames::uInterp, std::move(uInterp_));
    profile.set(ufo::ProfileVariableNames::vInterp, std::move(vInterp_));
    profile.set(ufo::ProfileVariableNames::LogP, std::move(LogP_));
    std::vector <int> NumStd(profile.getNumProfileLevels(), std::move(NumStd_));
    std::vector <int> NumSig(profile.getNumProfileLevels(), std::move(NumSig_));
    profile.set(ufo::ProfileVariableNames::NumStd, std::move(NumStd));
    profile.set(ufo::ProfileVariableNames::NumSig, std::move(NumSig));
  }
}  // namespace ufo
