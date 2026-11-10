/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include <set>
#include <utility>

#include "ufo/profile/ProfileCheckTime.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckTime>
  makerProfileCheckTime_("Time");

  ProfileCheckTime::ProfileCheckTime
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckTime::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Time check" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_surface_level_u,
       ufo::ProfileVariableNames::diagflags_perm_reject_u,
       ufo::ProfileVariableNames::diagflags_perm_reject_v};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileCheckTime::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsUSurfaceLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_u);
    std::vector <bool> &diagFlagsUPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_u);
    std::vector <bool> &diagFlagsVPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_v);

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::ObsType);
    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(ObsType,
                                         pressures,
                                         diagFlagsUSurfaceLevel,
                                         diagFlagsUPermReject,
                                         diagFlagsVPermReject)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Time checks will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(ObsType,
                                                    pressures,
                                                    diagFlagsUSurfaceLevel,
                                                    diagFlagsUPermReject,
                                                    diagFlagsVPermReject)
                         << std::endl;
      return;
    }

    // Reject sonde wind values for short period after launch.
    const float SondeLaunchWindRej = options_.TimeCheck_SondeLaunchWindRej.value();
    // Firstly determine surface pressure.
    float PSurf = 0.0;
    if (!diagFlagsUSurfaceLevel.empty() && SondeLaunchWindRej > 0.0 &&
        !ModelLevels &&
        (ObsType[0] != ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf)) {
      PSurf = pressures[0];
      for (size_t jlev = 0;
           jlev < std::min(static_cast<int>(numProfileLevels), 10);
           ++jlev) {
        if (diagFlagsUSurfaceLevel[jlev]) {
          PSurf = pressures[jlev];
          break;
        }
      }
    }

    // If surface pressure is nonzero, perform the wind rejection.
    if (PSurf > 0.0) {
      int NWindRej = 0;  // Number of wind levels rejected
      const float PLimit = PSurf - SondeLaunchWindRej * 100.0;  // Convert from hPa to Pa
      for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (pressures[jlev] > 0.0 && pressures[jlev] < PLimit) break;
        if (!diagFlagsUPermReject.empty()) diagFlagsUPermReject[jlev] = true;
        if (!diagFlagsVPermReject.empty()) diagFlagsVPermReject[jlev] = true;
        NWindRej++;
      }
      oops::Log::debug() << "Wind rejection: "
                         << "Psurf = " << PSurf * 0.01 << " hPa, "
                         << "NWindRej = " << NWindRej << std::endl;
    }
  }
}  // namespace ufo
