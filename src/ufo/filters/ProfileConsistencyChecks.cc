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

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/ProfileConsistencyCheckParameters.h"
#include "ufo/filters/ProfileConsistencyChecks.h"

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileChecker.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  // -----------------------------------------------------------------------------

  ProfileConsistencyChecks::ProfileConsistencyChecks
  (ioda::ObsSpace & obsdb,
   const eckit::Configuration & config,
   std::shared_ptr<ioda::ObsDataVector<int> > flags,
   std::shared_ptr<ioda::ObsDataVector<float> > obserr)
    : FilterBase(obsdb, config, flags, obserr)
  {
    options_.reset(new ProfileConsistencyCheckParameters());
    options_->deserialize(config);

    allvars_ += Variables(filtervars_, "HofX");

    // Throw exception if expected configuration option is missing.
    // It is essential for observations to be grouped according to (e.g.) station ID
    // (unless there is only one profile in the sample, which would be very unusual)
    if (obsdb.obs_group_var().empty())
      throw eckit::BadParameter("group variable is empty.", Here());
  }

  // -----------------------------------------------------------------------------

  ProfileConsistencyChecks::~ProfileConsistencyChecks() {}

  // -----------------------------------------------------------------------------

  void ProfileConsistencyChecks::applyFilter(const std::vector<bool> & apply,
                                             const Variables & filtervars,
                                             std::vector<std::vector<bool>> & flagged) const
  {
    print(oops::Log::trace());

    const int nlocs = static_cast <int> (obsdb_.nlocs());
    const int nprofs = static_cast <int> (obsdb_.nrecs());

    // Handles individual profile data
    ProfileDataHandler profileDataHandler(obsdb_,
                                          options_->DHParameters,
                                          apply);

    // (Optionally) validates check results against OPS values
    ProfileCheckValidator profileCheckValidator(*options_,
                                                profileDataHandler);

    // Applies checks to each profile
    ProfileChecker profileChecker(*options_,
                                  profileDataHandler,
                                  profileCheckValidator);

    // Loop over profiles
    oops::Log::debug() << "Starting loop over profiles..." << std::endl;

    for (int jprof = 0; jprof < nprofs; ++jprof) {
      oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

      /// Initialise the next profile prior to applying checks.
      profileDataHandler.initialiseNextProfile();

      // Print station ID if requested
      if (options_->PrintStationID.value()) {
        const std::vector <std::string> &station_ID =
          profileDataHandler.get<std::string>(ufo::VariableNames::station_ID);
        if (!station_ID.empty())
          oops::Log::debug() << "Station ID: " << station_ID[0] << std::endl;
      }

      // Run checks
      profileChecker.runChecks();

      // Update information, including the 'flagged' vector, for this profile.
      profileDataHandler.updateProfileInformation(filtervars.nvars(), flagged);

      // Optionally compare check results with OPS values
      if (options_->compareWithOPS.value() && profileChecker.getBasicCheckResult()) {
        profileCheckValidator.validate();
        nMismatches_.emplace_back(profileCheckValidator.getMismatches());
      }
    }

    // Write out any quantities that may have changed to obsdb
    profileDataHandler.writeQuantitiesToObsdb();

    oops::Log::debug() << "... Finished loop over profiles" << std::endl;
    oops::Log::debug() << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ProfileConsistencyChecks::print(std::ostream & os) const {
    os << "ProfileConsistencyChecks: config = " << config_ << std::endl;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo
