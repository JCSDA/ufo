/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
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

#include "ufo/profile/ProfileChecker.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileData.h"
#include "ufo/profile/ProfileFlags.h"
#include "ufo/profile/ProfileIndices.h"

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

    // Throw exception if expected configuration options are missing.
    // It is essential for observations to be grouped according to (e.g.) station ID
    // (unless there is only one profile in the sample, which would be very unusual)
    if (obsdb.obs_group_var().empty())
      throw eckit::BadParameter("obsdatain.obsgrouping.group_var is empty.", Here());
  }

  // -----------------------------------------------------------------------------

  ProfileConsistencyChecks::~ProfileConsistencyChecks() {}

  // -----------------------------------------------------------------------------

  void ProfileConsistencyChecks::applyFilter(const std::vector<bool> & apply,
                                             const Variables & filtervars,
                                             std::vector<std::vector<bool>> & flagged) const
  {
    const int nlocs = static_cast <int> (obsdb_.nlocs());
    const int nprofs = static_cast <int> (obsdb_.nrecs());

    // Determines indices of profile's observations in entire sample
    ProfileIndices profileIndices(obsdb_,
                                  *options_,
                                  apply);

    // Gets individual profile data
    ProfileData profileData(obsdb_,
                            *options_,
                            profileIndices);

    // Gets individual profile flags and counters and modifies them as required
    ProfileFlags profileFlags(obsdb_,
                              *options_,
                              profileIndices);

    // (Optionally) validates check results against OPS values
    ProfileCheckValidator profileCheckValidator(obsdb_,
                                                *options_,
                                                profileIndices);

    // Applies checks to each profile
    const ProfileChecker profileChecker(*options_,
                                        profileIndices,
                                        profileData,
                                        profileFlags,
                                        profileCheckValidator);

    // Loop over profiles
    oops::Log::debug() << "Starting loop over profiles..." << std::endl;

    for (int jprof = 0; jprof < nprofs; ++jprof) {
      oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

      // Determine indices in total sample that correspond to this profile
      profileIndices.determineProfileIndices();

      // Load values of physical variables (for both observations and model)
      profileData.fillProfileValues();
      oops::Log::debug() << "Station ID: " << profileData.getStationID() << std::endl;

      // Load QC flags and counters
      profileFlags.fillProfileValues();
      profileFlags.setProfileNum(jprof);

      // Run checks
      profileChecker.runChecks();

      // Set final report flags if a certain number of errors have occurred
      profileFlags.setFinalReportFlags();

      // Update flags in the entire sample (using values in this profile)
      profileFlags.updateFlags();

      // Optionally compare check results with OPS values
      if (options_->compareWithOPS.value() && profileFlags.getBasicCheckResult()) {
        profileCheckValidator.setReportFlags(profileFlags.getReportFlags());
        profileCheckValidator.fillProfileValues();
        profileCheckValidator.setProfileNum(jprof);
        profileCheckValidator.validate();
        nMismatches_.emplace_back(profileCheckValidator.getMismatches());
      }
    }

    // Modify flagged
    profileFlags.setFlagged(nlocs, filtervars.nvars(), flagged);

    // Write out flags, counters and data corrections to obsdb
    profileFlags.writeFlags();

    oops::Log::debug() << "... Finished loop over profiles" << std::endl;
    oops::Log::debug() << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ProfileConsistencyChecks::print(std::ostream & os) const {
    os << "ProfileConsistencyChecks: config = " << config_ << std::endl;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo
