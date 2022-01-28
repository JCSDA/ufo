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
#include <utility>
#include <vector>

#include "ufo/filters/ConventionalProfileProcessingParameters.h"

namespace ufo {
  class ProfileDataHandler;
}

namespace ufo {

  /// Information on each subgroup of checks.
  struct CheckSubgroup {
    /// \p runOnEntireSample specifies whether the checks in this subgroup run on all
    /// profiles at once.
    bool runOnEntireSample;
    /// \p checkNames contains the names of the checks in this subgroup.
    std::vector<std::string> checkNames;
  };

  /// \brief Profile QC checker
  ///
  /// Runs the various QC checks on individual profiles and modifies flags accordingly.
  class ProfileChecker {
   public:
    explicit ProfileChecker(const ConventionalProfileProcessingParameters &options);

    /// Type for container of check subgroups.
    typedef std::vector <CheckSubgroup> CheckSubgroupList;

    /// Run all checks requested
    void runChecks(ProfileDataHandler &profileDataHandler,
                   const CheckSubgroup &subGroupChecks);

    /// Get basic check result
    bool getBasicCheckResult() {return basicCheckResult_;}

    /// Set basic check result
    void setBasicCheckResult(bool result) {basicCheckResult_ = result;}

    /// Get container of check subgroups
    CheckSubgroupList getCheckSubgroups() {return checkSubgroups_;}

    /// Indicate whether at least one check requires HofX to have been calculated.
    bool requiresHofX() const {return requiresHofX_;}

    /// Get vector of GeoVaL names for all checks.
    oops::Variables getGeoVaLNames() const {return GeoVaLNames_;}

    /// Get vector of validation GeoVaL names for all checks.
    oops::Variables getValidationGeoVaLNames() const {return validationGeoVaLNames_;}

    /// Get vector of obs diagnostic names for all checks.
    oops::Variables getObsDiagNames() const {return obsDiagNames_;}

   private:
    /// Configurable parameters
    const ConventionalProfileProcessingParameters &options_;

    /// Checks to perform
    std::vector <std::string> checks_;

    /// Basic check result
    bool basicCheckResult_ = true;

    /// Subgroups of checks with the same mode of operation.
    CheckSubgroupList checkSubgroups_;

    /// Variable indicating whether at least one check requires HofX to have been calculated.
    bool requiresHofX_;

    /// Names of all required GeoVaLs.
    oops::Variables GeoVaLNames_;

    /// Names of all validation GeoVaLs.
    oops::Variables validationGeoVaLNames_;

    /// Names of all required obs diagnostics.
    oops::Variables obsDiagNames_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKER_H_
