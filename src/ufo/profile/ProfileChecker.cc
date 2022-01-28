/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <set>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckBasic.h"
#include "ufo/profile/ProfileChecker.h"
#include "ufo/profile/ProfileCheckHydrostatic.h"
#include "ufo/profile/ProfileCheckInterpolation.h"
#include "ufo/profile/ProfileCheckRH.h"
#include "ufo/profile/ProfileCheckSamePDiffT.h"
#include "ufo/profile/ProfileCheckSign.h"
#include "ufo/profile/ProfileCheckUInterp.h"
#include "ufo/profile/ProfileCheckUnstableLayer.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {
  ProfileChecker::ProfileChecker(const ConventionalProfileProcessingParameters &options)
    : options_(options),
      checks_(options.Checks.value()),
      requiresHofX_(false)
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

    // Produce check subgroups.
    // A subgroup is the longest sequence of consecutive checks which have
    // the same mode of operation.
    bool isFirst = true;  // First check considered; used to initialise checkMode.
    bool checkMode = true;  // Mode of operation of the current check.
    std::vector <std::string> checkNames;  // Filled anew for each subgroup.
    // Also fill a set of any required GeoVaL names.
    for (const auto& check : checks_) {
      // Instantiate each check and check its mode of operation.
      std::unique_ptr<ProfileCheckBase> profileCheck =
        ProfileCheckFactory::create(check,
                                    options_);
      if (profileCheck) {
        if (isFirst) {
          checkMode = profileCheck->runOnEntireSample();
          isFirst = false;
        }
        if (profileCheck->runOnEntireSample() == checkMode) {
          checkNames.push_back(check);
        } else {
          checkSubgroups_.push_back({checkMode, checkNames});
          checkNames.clear();
          checkNames.push_back(check);
          // Invert checkMode whenever a check with a different mode is reached.
          checkMode = !checkMode;
        }
        requiresHofX_ = requiresHofX_ || profileCheck->requiresHofX();
        GeoVaLNames_ += profileCheck->getGeoVaLNames();
        validationGeoVaLNames_ += profileCheck->getValidationGeoVaLNames();
        obsDiagNames_ += profileCheck->getObsDiagNames();
      } else {
        throw eckit::NotImplemented("Have not implemented a check for " + check, Here());
      }
    }
    // Fill checkSubgroups with the final list to be produced.
    checkSubgroups_.push_back({checkMode, checkNames});
  }

  void ProfileChecker::runChecks(ProfileDataHandler &profileDataHandler,
                                 const CheckSubgroup &subGroupChecks)
  {
    // Run all checks requested
    for (const auto& check : subGroupChecks.checkNames) {
      std::unique_ptr<ProfileCheckBase> profileCheck =
        ProfileCheckFactory::create(check,
                                    options_);
      if (profileCheck) {
        // Ensure correct type of check has been requested.
        if (profileCheck->runOnEntireSample() == subGroupChecks.runOnEntireSample) {
          // For checks on the entire sample, reset profile indices
          // prior to looping through the profiles
          if (profileCheck->runOnEntireSample())
            profileDataHandler.resetProfileIndices();
          // Run check
          profileCheck->runCheck(profileDataHandler);
          // Actions taken if a single profile was processed.
          if (!profileCheck->runOnEntireSample()) {
            // Fill validation information if required
            if (options_.compareWithOPS.value())
              profileCheck->fillValidationData(profileDataHandler);
            // Do not proceed if basic checks failed
            if (!profileCheck->getResult() && check == "Basic") {
              oops::Log::debug() << "Basic checks failed" << std::endl;
              setBasicCheckResult(false);
              break;
            }
          }
        }
      } else {
        throw eckit::NotImplemented("Have not implemented a check for " + check, Here());
      }
    }
  }
}  // namespace ufo
