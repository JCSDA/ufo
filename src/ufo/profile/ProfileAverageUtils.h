/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGEUTILS_H_
#define UFO_PROFILE_PROFILEAVERAGEUTILS_H_

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oops/util/missingValues.h"

#include "ufo/profile/ProfileDataHolder.h"

namespace ufo {
  class ProfileDataHandler;
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

    /// Set values in a profile to missing.
    template <typename T>
    static void setProfileMissing
      (ProfileDataHolder & profile,
       const std::vector <std::string> & variableNames)
    {
      const T missing = util::missingValue(missing);
      for (const std::string & variableName : variableNames) {
        profile.set<T>(variableName,
                       std::move(std::vector<T>(profile.getNumProfileLevels(), missing)));
      }
    }

    /// Transfer one variable in a profile to another variable in the same profile.
    template <typename T>
      static void copyProfileValues(ProfileDataHolder & profile,
                                    const std::string & sourceVariableName,
                                    const std::string & destinationVariableName)
    {
      profile.set<T>(destinationVariableName,
                     std::move(profile.get<T>(sourceVariableName)));
    }
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEUTILS_H_
