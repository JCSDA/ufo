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

namespace ufo {
  class ProfileCheckValidator;
  class ProfileDataHandler;
}

namespace ufo {

  /// \brief Profile QC checker
  ///
  /// Runs the various QC checks on individual profiles and modifies flags accordingly.
  class ProfileChecker {
   public:
    ProfileChecker(const ProfileConsistencyCheckParameters &options,
                   ProfileDataHandler &profileDataHandler,
                   ProfileCheckValidator &profileCheckValidator);

    /// Run all checks requested
    void runChecks();

    /// Get basic check result
    bool getBasicCheckResult() {return basicCheckResult_;}

    /// Set basic check result
    void setBasicCheckResult(bool result) {basicCheckResult_ = result;}

   private:
    /// Configurable parameters
    const ProfileConsistencyCheckParameters &options_;

    /// Profile data
    ProfileDataHandler &profileDataHandler_;

    /// Profile check validator
    ProfileCheckValidator &profileCheckValidator_;

    /// Checks to perform
    std::vector <std::string> checks_;

    /// Basic check result
    bool basicCheckResult_ = true;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKER_H_
