/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckHistory.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckHistory>
  makerProfileCheckHistory_("History");

  ProfileCheckHistory::ProfileCheckHistory
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckHistory::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " History check" << std::endl;

    const size_t nprofs = profileDataHandler.getObsdb().nrecs();

    // Vector of profiles containing data for the history check.
    std::vector <ProfileDataHistory> profiles;
    for (size_t jprof = 0; jprof < nprofs; ++jprof) {
      profileDataHandler.initialiseNextProfile();
      ProfileDataHistory profile(profileDataHandler);
      profile.fill();
      profiles.emplace_back(profile);
    }

    // todo: this is how to access data within individual profiles
    /*
    for (auto &profile : profiles) {
      const std::vector<float> &U = profile.get<float>(ufo::VariableNames::obs_eastward_wind);
    }
    */
  }
}  // namespace ufo
