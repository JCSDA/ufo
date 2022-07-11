/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ProfileLevelCount.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<ProfileLevelCount>
                       makerProfileLevelCount_("ProfileLevelCount");

// -----------------------------------------------------------------------------

ProfileLevelCount::ProfileLevelCount(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  invars_ += getAllWhereVariables(options_.where);
}

// -----------------------------------------------------------------------------

ProfileLevelCount::~ProfileLevelCount() {}

// -----------------------------------------------------------------------------

void ProfileLevelCount::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "ProfileLevelCount::compute started" << std::endl;

  // ObsSpace.
  const ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into profiles.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Number of locations.
  const size_t nlocs = obsdb.nlocs();

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<std::size_t> &recnums = obsdb.recidx_all_recnums();

  // Number of profiles in the ObsSpace.
  const std::size_t nprofs = recnums.size();

  // Vector of locations that pass the 'where' clause in the sample
  // (all true if there is no where clause).
  const std::vector<bool> apply = processWhere(options_.where, in, options_.whereOperator);

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

    // Get locations of this profile.
    const std::vector<std::size_t> &locs = obsdb.recidx_vector(recnums[jprof]);

    // Count number of valid locations in this profile.
    int count = 0;
    for (std::size_t loc : locs)
      if (apply[loc])
        count++;

    // Write out count to all locations in output vector.
    for (std::size_t loc : locs)
      out[0][loc] = count;
  }

  oops::Log::trace() << "ProfileLevelCount::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ProfileLevelCount::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
