/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckUInterpAlternative.h"
#include "ufo/profile/VariableNames.h"

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
      {ufo::VariableNames::qcflags_eastward_wind,
       ufo::VariableNames::counter_NumSamePErrObs,
       ufo::VariableNames::counter_NumInterpErrObs,
       ufo::VariableNames::extended_obs_space};
    std::vector <std::string> variableNamesFloat =
      {ufo::VariableNames::obs_air_pressure,
       ufo::VariableNames::obs_eastward_wind,
       ufo::VariableNames::obs_northward_wind};
    if (options_.compareWithOPS.value()) {
      variableNamesInt.insert(variableNamesInt.end(),
                              {ufo::VariableNames::StdLev,
                                  ufo::VariableNames::SigAbove,
                                  ufo::VariableNames::SigBelow,
                                  ufo::VariableNames::LevErrors,
                                  ufo::VariableNames::NumStd,
                                  ufo::VariableNames::NumSig});
      variableNamesFloat.insert(variableNamesFloat.end(),
                                {ufo::VariableNames::uInterp,
                                    ufo::VariableNames::vInterp,
                                    ufo::VariableNames::LogP});
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
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
      const auto &extended_obs_space = profile.get<int>(ufo::VariableNames::extended_obs_space);
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
      profile.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &uObs =
      profile.get<float>(ufo::VariableNames::obs_eastward_wind);
    const std::vector <float> &vObs =
      profile.get<float>(ufo::VariableNames::obs_northward_wind);
    std::vector <int> &uFlags =
      profile.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &NumSamePErrObs =
      profile.get<int>(ufo::VariableNames::counter_NumSamePErrObs);
    std::vector <int> &NumInterpErrObs =
      profile.get<int>(ufo::VariableNames::counter_NumInterpErrObs);

    if (pressures.empty() ||
        pressures.front() > options_.BChecks_maxValidP.value() ||
        pressures.back() < options_.BChecks_minValidP.value()) {
      return;
    }

    if (!oops::allVectorsSameNonZeroSize(pressures, uObs, vObs, uFlags)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, uObs, vObs, uFlags)
                         << std::endl;
      return;
    }

    // Populate the standard level diagnostic flags from the OPS-style QC flags.
    // todo(ctgh): modify this when the ProfileDataHolder class can accept a vector of bools.
    std::vector <bool> uDiagFlagsProfileStdLev;
    std::transform(uFlags.begin(), uFlags.end(),
                   std::back_inserter(uDiagFlagsProfileStdLev),
                   [](int flag){return flag & ufo::MetOfficeQCFlags::Profile::StandardLevelFlag;});

    calcStdLevelsUV(numProfileLevels, pressures, uObs, vObs, uDiagFlagsProfileStdLev);

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

  void ProfileCheckUInterpAlternative::fillValidationData(ProfileDataHolder &profile)
  {
    profile.set(ufo::VariableNames::StdLev, std::move(StdLev_));
    profile.set(ufo::VariableNames::SigAbove, std::move(SigAbove_));
    profile.set(ufo::VariableNames::SigBelow, std::move(SigBelow_));
    profile.set(ufo::VariableNames::LevErrors, std::move(LevErrors_));
    profile.set(ufo::VariableNames::uInterp, std::move(uInterp_));
    profile.set(ufo::VariableNames::vInterp, std::move(vInterp_));
    profile.set(ufo::VariableNames::LogP, std::move(LogP_));
    std::vector <int> NumStd(profile.getNumProfileLevels(), std::move(NumStd_));
    std::vector <int> NumSig(profile.getNumProfileLevels(), std::move(NumSig_));
    profile.set(ufo::VariableNames::NumStd, std::move(NumStd));
    profile.set(ufo::VariableNames::NumSig, std::move(NumSig));
  }
}  // namespace ufo
