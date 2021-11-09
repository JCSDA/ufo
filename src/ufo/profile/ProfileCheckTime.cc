/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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

  void ProfileCheckTime::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Time check" << std::endl;

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::VariableNames::ObsType);
    const std::vector <float> &level_time =
       profileDataHandler.get<float>(ufo::VariableNames::obs_level_time);
    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    std::vector <int> &uFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_northward_wind);
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::VariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(ObsType, pressures, uFlags, vFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Time checks will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(ObsType, pressures, uFlags, vFlags)
                           << std::endl;
      return;
    }

    // Flag any observations that appear outside of the current time window.
    // The variable level_time is equal to the number of seconds relative to
    // the middle of the time window. level_time is compared to half of the
    // assimilation window length.
    const float halfWindowLength = 0.5 * (profileDataHandler.getObsdb().windowEnd() -
                                          profileDataHandler.getObsdb().windowStart()).toSeconds();
    // Time QC flags. A value of true indicates an observation
    // lies outside the time window.
    // These flags are used in subsequent background checks.
    std::vector <int> timeFlags(numProfileLevels, false);
    if (!level_time.empty() && !ModelLevels) {
      for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
        const float leveltime = level_time[jlev];
        if (leveltime == missingValueFloat) continue;
        timeFlags[jlev] = (leveltime < (-halfWindowLength - 0.5) ||
                           leveltime > (halfWindowLength - 0.5));
      }
    }

    // Reject sonde wind values for short period after launch.
    const float SondeLaunchWindRej = options_.TimeCheck_SondeLaunchWindRej.value();
    // Firstly determine surface pressure.
    float PSurf = 0.0;
    if (!uFlags.empty() && SondeLaunchWindRej > 0.0 &&
        !ModelLevels &&
        (ObsType[0] != ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf)) {
      PSurf = pressures[0];
      for (size_t jlev = 0;
           jlev < std::min(static_cast<int>(numProfileLevels), 10);
           ++jlev) {
        if (uFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag) {
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
        if (!uFlags.empty()) uFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
        if (!vFlags.empty()) vFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
        NWindRej++;
      }
      oops::Log::debug() << "Wind rejection: "
                         << "Psurf = " << PSurf * 0.01 << " hPa, "
                         << "NWindRej = " << NWindRej << std::endl;
    }

    // Store the time flags for use in later checks.
    profileDataHandler.set<int>(ufo::VariableNames::qcflags_time, std::move(timeFlags));
  }
}  // namespace ufo
