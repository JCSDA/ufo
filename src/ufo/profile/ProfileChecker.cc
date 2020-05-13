/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckBasic.h"
#include "ufo/profile/ProfileChecker.h"
#include "ufo/profile/ProfileCheckHydrostatic.h"
#include "ufo/profile/ProfileCheckInterpolation.h"
#include "ufo/profile/ProfileCheckSamePDiffT.h"
#include "ufo/profile/ProfileCheckSign.h"
#include "ufo/profile/ProfileCheckUnstableLayer.h"

#include "ufo/utils/Flags.h"

namespace ufo {
  ProfileChecker::ProfileChecker(const ProfileConsistencyCheckParameters &options,
                                 const ProfileIndices &profileIndices,
                                 const ProfileData &profileData,
                                 ProfileFlags &profileFlags,
                                 ProfileCheckValidator &profileCheckValidator)
    : options_(options),
      profileIndices_(profileIndices),
      profileData_(profileData),
      profileFlags_(profileFlags),
      profileCheckValidator_(profileCheckValidator),
      checks_(options.Checks.value())
  {
    // Ensure basic checks are always performed first
    auto it_checks = std::find(checks_.begin(), checks_.end(), "Basic");
    if (it_checks == checks_.end()) {
      // If "Basic" not present, insert at the start
      checks_.insert(checks_.begin(), "Basic");
    } else if (it_checks != checks_.begin()) {
      // If "Basic" is present but not at the start, move it there
      std::rotate(checks_.begin(), it_checks, it_checks + 1);
    }
  }

  void ProfileChecker::runChecks() const
  {
    // Run all checks requested
    for (auto check : checks_) {
      std::unique_ptr <ProfileCheckBase> profileCheck(nullptr);
      if (check == "Basic") {
        profileCheck.reset(new ProfileCheckBasic(options_,
                                                 profileIndices_,
                                                 profileData_,
                                                 profileFlags_,
                                                 profileCheckValidator_));
      } else if (check == "SamePDiffT") {
        profileCheck.reset(new ProfileCheckSamePDiffT(options_,
                                                      profileIndices_,
                                                      profileData_,
                                                      profileFlags_,
                                                      profileCheckValidator_));
      } else if (check == "Sign") {
        profileCheck.reset(new ProfileCheckSign(options_,
                                                profileIndices_,
                                                profileData_,
                                                profileFlags_,
                                                profileCheckValidator_));
      } else if (check == "UnstableLayer") {
        profileCheck.reset(new ProfileCheckUnstableLayer(options_,
                                                         profileIndices_,
                                                         profileData_,
                                                         profileFlags_,
                                                         profileCheckValidator_));
      } else if (check == "Interpolation") {
        profileCheck.reset(new ProfileCheckInterpolation(options_,
                                                         profileIndices_,
                                                         profileData_,
                                                         profileFlags_,
                                                         profileCheckValidator_));
      } else if (check == "Hydrostatic") {
        profileCheck.reset(new ProfileCheckHydrostatic(options_,
                                                       profileIndices_,
                                                       profileData_,
                                                       profileFlags_,
                                                       profileCheckValidator_));
      }
      if (profileCheck) {
        profileCheck->runCheck();
        // Fill validation information if required
        if (options_.compareWithOPS.value()) {
          profileCheck->fillValidator();
        }
        // Do not proceed if basic checks failed
        if (check == "Basic" && !profileCheck->getResult()) {
          oops::Log::debug() << "Basic checks failed" << std::endl;
          profileFlags_.setBasicCheckResult(false);
          break;
        }
      } else {
        throw eckit::NotImplemented("Have not implemented a check for " + check, Here());
      }
    }
  }
}  // namespace ufo
