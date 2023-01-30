/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/ConventionalProfileProcessing.h"
#include "ufo/filters/ConventionalProfileProcessingParameters.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

  // -----------------------------------------------------------------------------

  ConventionalProfileProcessing::ConventionalProfileProcessing
  (ioda::ObsSpace & obsdb,
   const Parameters_ & parameters,
   std::shared_ptr<ioda::ObsDataVector<int> > flags,
   std::shared_ptr<ioda::ObsDataVector<float> > obserr)
    : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
  {
    // Determine whether any check requires HofX to have been calculated.
    // If so, add the filter variables to \c allvars_.
    // This is performed here because \c allvars_ must be filled prior to the filter being run.
    // The \c profileChecker class must be instantiated in order to do this.
    const ProfileChecker profileChecker(options_);
    if (profileChecker.requiresHofX())
      allvars_ += Variables(filtervars_, "HofX");

    // Add the names of any required GeoVaLs to \c allvars_.
    // The names of the required GeoVaLs are listed in the \c getGeoVaLNames() function
    // of each check (which returns an empty set by default).
    const auto &requiredGeoVaLNames = profileChecker.getGeoVaLNames();
    if (requiredGeoVaLNames.size() > 0) {
      ufo::Variables reqGeoVaLs(requiredGeoVaLNames);
      allvars_ += Variables(reqGeoVaLs, "GeoVaLs");
    }

    // Add the names of any required obs diagnostics to \c allvars_.
    // The names of the required obs diagnostics are listed in the \c getObsDiagNames() function
    // of each check (which returns an empty set by default).
    const auto &requiredObsDiagNames = profileChecker.getObsDiagNames();
    if (requiredObsDiagNames.size() > 0) {
      ufo::Variables reqObsDiags(requiredObsDiagNames);
      allvars_ += Variables(reqObsDiags, "ObsDiag");
    }

    if (options_.compareWithOPS.value()) {
      const auto &validationGeoVaLNames = profileChecker.getValidationGeoVaLNames();
      if (validationGeoVaLNames.size() > 0) {
        ufo::Variables valGeoVaLs(validationGeoVaLNames);
        allvars_ += Variables(valGeoVaLs, "GeoVaLs");
      }
    }

    // It is essential for observations to be grouped according to (e.g.) station ID
    // (unless there is only one profile in the sample, which would be very unusual).
    // Throw an exception if the "group variable" configuration option is missing.
    if (obsdb.obs_group_vars().empty())
      throw eckit::BadParameter("group variables configuration is empty.", Here());
  }

  // -----------------------------------------------------------------------------

  ConventionalProfileProcessing::~ConventionalProfileProcessing() {}

  // -----------------------------------------------------------------------------

  void ConventionalProfileProcessing::individualProfileChecks
  (ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator,
   ProfileChecker &profileChecker,
   const CheckSubgroup &subGroupChecks) const
  {
    const int nprofs = static_cast <int> (obsdb_.nrecs());

    // Reset profile indices prior to looping through entire sample.
    profileDataHandler.resetProfileIndices();

    // Loop over profiles
    oops::Log::debug() << "Starting loop over profiles..." << std::endl;

    for (int jprof = 0; jprof < nprofs; ++jprof) {
      oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

      // Initialise the next profile prior to applying checks.
      profileDataHandler.initialiseNextProfile();

      // Print station ID if requested
      if (options_.PrintStationID.value()) {
        const std::vector <std::string> &station_ID =
          profileDataHandler.get<std::string>(ufo::VariableNames::station_ID);
        if (!station_ID.empty())
          oops::Log::debug() << "Station ID: " << station_ID[0] << std::endl;
      }

      // Run checks
      profileChecker.runChecks(profileDataHandler,
                               subGroupChecks);

      // Update information, including the 'flagged' vector, for this profile.
      profileDataHandler.updateProfileInformation();
    }

    // Write various quantities to the obsdb.
    profileDataHandler.writeQuantitiesToObsdb();

    oops::Log::debug() << "... Finished loop over profiles" << std::endl;
    oops::Log::debug() << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ConventionalProfileProcessing::entireSampleChecks
  (ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator,
   ProfileChecker &profileChecker,
   const CheckSubgroup &subGroupChecks) const
  {
    oops::Log::debug() << "Running checks on entire profile sample..." << std::endl;

    // Run checks
    profileChecker.runChecks(profileDataHandler,
                             subGroupChecks);

    // Write various quantities to the obsdb.
    profileDataHandler.writeQuantitiesToObsdb();

    oops::Log::debug() << "... Finished running checks" << std::endl;
    oops::Log::debug() << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ConventionalProfileProcessing::applyFilter
  (const std::vector<bool> & apply,
   const Variables & filtervars,
   std::vector<std::vector<bool>> & flagged) const
  {
    print(oops::Log::trace());

    // Handles individual profile data
    ProfileDataHandler profileDataHandler(data_,
                                          *flags_,
                                          options_.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);

    // (Optionally) validates check results against OPS values
    ProfileCheckValidator profileCheckValidator(options_);

    // Applies checks to each profile
    ProfileChecker profileChecker(options_);

    // Loop over each check subgroup in turn.
    const auto checkSubgroups = profileChecker.getCheckSubgroups();
    for (const auto& checkSubgroup : checkSubgroups) {
      if (checkSubgroup.runOnEntireSample) {
        // Run checks that use all of the profiles at once.
        entireSampleChecks(profileDataHandler,
                           profileCheckValidator,
                           profileChecker,
                           checkSubgroup);
      } else {
        // Run checks on individual profiles sequentially.
        individualProfileChecks(profileDataHandler,
                                profileCheckValidator,
                                profileChecker,
                                checkSubgroup);
      }
    }

    // Optionally compare check results with OPS values
    if (options_.compareWithOPS.value()) {
      oops::Log::debug() << " Comparing values against OPS equivalents..." << std::endl;
      profileDataHandler.resetProfileIndices();
      for (size_t jprof = 0; jprof < obsdb_.nrecs(); ++jprof) {
        oops::Log::debug() << " Profile " << jprof + 1 << std::endl;
        profileDataHandler.initialiseNextProfile();
        profileCheckValidator.validate(profileDataHandler, obsdb_.comm().size());
        nMismatches_.emplace_back(profileCheckValidator.getMismatches());
        oops::Log::debug() << std::endl;
      }
    }
  }

  // -----------------------------------------------------------------------------

  void ConventionalProfileProcessing::print(std::ostream & os) const {
    os << "ConventionalProfileProcessing: config = " << options_ << std::endl;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo
