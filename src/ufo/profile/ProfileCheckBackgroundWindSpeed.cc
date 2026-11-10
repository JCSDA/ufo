/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckBackgroundWindSpeed.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundWindSpeed>
  makerProfileCheckBackgroundWindSpeed_("BackgroundWindSpeed");

  ProfileCheckBackgroundWindSpeed::ProfileCheckBackgroundWindSpeed
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBackgroundWindSpeed::doPrepareOutputColumns(
      ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Background check for wind velocity" << std::endl;

    const std::vector<std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::pge_eastward_wind,
       ufo::ProfileVariableNames::pge_northward_wind};
    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_interpolation_u,
       ufo::ProfileVariableNames::diagflags_back_perf_u,
       ufo::ProfileVariableNames::diagflags_back_reject_u,
       ufo::ProfileVariableNames::diagflags_perm_reject_u,
       ufo::ProfileVariableNames::diagflags_final_reject_u,
       ufo::ProfileVariableNames::diagflags_interpolation_v,
       ufo::ProfileVariableNames::diagflags_back_perf_v,
       ufo::ProfileVariableNames::diagflags_back_reject_v,
       ufo::ProfileVariableNames::diagflags_perm_reject_v,
       ufo::ProfileVariableNames::diagflags_final_reject_v};

    for (const auto &variableFloat : variableNamesFloat) {
      profileDataHandler.get<float>(variableFloat);
    }
    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileCheckBackgroundWindSpeed::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <float> &uPGE =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::pge_eastward_wind);
    std::vector <bool> &diagFlagsUInterp =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_u);
    std::vector <bool> &diagFlagsUBackPerf =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_perf_u);
    std::vector <bool> &diagFlagsUBackReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_reject_u);
    std::vector <bool> &diagFlagsUPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_u);
    std::vector <bool> &diagFlagsUFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_u);
    std::vector <float> &vPGE =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::pge_northward_wind);
    std::vector <bool> &diagFlagsVInterp =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_v);
    std::vector <bool> &diagFlagsVBackPerf =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_perf_v);
    std::vector <bool> &diagFlagsVBackReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_back_reject_v);
    std::vector <bool> &diagFlagsVPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_v);
    std::vector <bool> &diagFlagsVFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_v);
    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &uObs =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_eastward_wind);
    const std::vector <float> &uObsErr =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obserr_eastward_wind);
    const std::vector <float> &uBkg =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_eastward_wind);
    const std::vector <float> &uBkgErr =
      profileDataHandler.get<float>
      (options_.bkgErrGroup.value() + "/" + options_.bkgErrName_eastward_wind.value());
    const std::vector <float> &vObs =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_northward_wind);
    const std::vector <float> &vObsErr =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obserr_northward_wind);
    const std::vector <float> &vBkg =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_northward_wind);
    const std::vector <float> &vBkgErr =
      profileDataHandler.get<float>
      (options_.bkgErrGroup.value() + "/" + options_.bkgErrName_northward_wind.value());
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();


    if (!oops::allVectorsSameNonZeroSize(uObs, uObsErr, uBkg, uBkgErr,
                                         uPGE,
                                         diagFlagsUInterp,
                                         vObs, vObsErr, vBkg, vBkgErr,
                                         vPGE,
                                         diagFlagsVInterp)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(uObs, uObsErr, uBkg, uBkgErr,
                                                    uPGE, diagFlagsUInterp,
                                                    vObs, vObsErr, vBkg, vBkgErr,
                                                    vPGE, diagFlagsVInterp)
                         << std::endl;
      return;
    }

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, options_.BkCheck_PdBad_uv.value());

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (diagFlagsUInterp[jlev])
        uPGE[jlev] = 0.5 + 0.5 * uPGE[jlev];
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           uObs,
                           uObsErr,
                           uBkg,
                           uBkgErr,
                           PdBad,
                           ModelLevels,
                           diagFlagsUBackPerf,
                           diagFlagsUBackReject,
                           diagFlagsUPermReject,
                           diagFlagsUFinalReject,
                           uPGE,
                           -1,
                           &vObs,
                           &vBkg);

    // Update v PGE and flags.
    vPGE = uPGE;
    diagFlagsVInterp = diagFlagsUInterp;
    diagFlagsVBackPerf = diagFlagsUBackPerf;
    diagFlagsVBackReject = diagFlagsUBackReject;
    diagFlagsVPermReject = diagFlagsUPermReject;
    diagFlagsVFinalReject = diagFlagsUFinalReject;
  }
}  // namespace ufo
