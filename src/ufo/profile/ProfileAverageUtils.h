/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGEUTILS_H_
#define UFO_PROFILE_PROFILEAVERAGEUTILS_H_

#include <string>
#include <vector>

namespace ufo {
  class ProfileDataHandler;
  class ProfileDataHolder;
}

namespace ufo {

  class ProfileAverageUtils {
   public:
    /// Modify filter QC flags based on values of the averaged data.
    static void passNonMissingAveragedObservations
      (ProfileDataHandler & profileDataHandler,
       std::vector <ProfileDataHolder> & profiles,
       const std::string & flag_name,
       const std::string & obs_name);

    /// Fill validation data for use in comparisons with OPS.
    static void fillValidationData
      (ProfileDataHolder & profile,
       bool extended_obs_space,
       const std::string & average_name,
       const std::string & qcflags_name,
       const std::string & geovals_testreference_name,
       const std::string & geovals_qcflags_name);
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEUTILS_H_
