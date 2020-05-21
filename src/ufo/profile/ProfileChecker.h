/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKER_H_
#define UFO_PROFILE_PROFILECHECKER_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/ProfileFlags.h"

namespace ufo {
  class ProfileCheckValidator;
  class ProfileData;
  class ProfileIndices;
}

namespace ufo {

  /// \brief Profile QC checker
  ///
  /// Runs the various QC checks on individual profiles and modifies flags accordingly.
  class ProfileChecker {
   public:
    ProfileChecker(const ProfileConsistencyCheckParameters &options,
                   const ProfileIndices &profileIndices,
                   const ProfileData &profileData,
                   ProfileFlags &profileFlags,
                   ProfileCheckValidator &profileCheckValidator);

    /// Run all checks requested
    void runChecks() const;

   private:
    /// Configurable parameters
    const ProfileConsistencyCheckParameters &options_;

    /// Indices of profile's observations in the entire sample
    const ProfileIndices &profileIndices_;

    /// Profile data
    const ProfileData &profileData_;

    /// Profile flags
    ProfileFlags &profileFlags_;

    /// Profile check validator
    ProfileCheckValidator &profileCheckValidator_;

    /// Checks to perform
    std::vector <std::string> checks_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKER_H_
