/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/FlagOriginalAndAveragedProfiles.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<FlagOriginalAndAveragedProfiles>

makerFlagOriginalAndAveragedProfiles_("flag original and averaged profiles");

// -----------------------------------------------------------------------------

FlagOriginalAndAveragedProfiles::FlagOriginalAndAveragedProfiles
(const FlagOriginalAndAveragedProfilesParameters &parameters)
  : allvars_(), parameters_(parameters) {
}

// -----------------------------------------------------------------------------

void FlagOriginalAndAveragedProfiles::apply(const Variables & vars,
                                            const std::vector<std::vector<bool>> & flagged,
                                            ObsFilterData & data,
                                            int filterQCflag,
                                            ioda::ObsDataVector<int> & flags,
                                            ioda::ObsDataVector<float> &) const {
  // ObsSpace.
  const ioda::ObsSpace & obsdb = data.obsspace();

  // Ensure observations have been grouped into profiles.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb.has("MetaData", "extendedObsSpace"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<std::size_t> &recnums = data.obsspace().recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Loop over each filter variable
  for (size_t jv = 0; jv < vars.nvars(); ++jv) {
    const size_t iv = flags.varnames().find(vars.variable(jv).variable());
    // Loop over original profiles.
    for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
      // Get locations of profile in the original ObsSpace and
      // the corresponding profile in the extended ObsSpace.
      // Assuming the extended ObsSpace has been configured correctly, which is
      // checked above, the profile in the extended ObsSpace is always located
      // nprofs positions further on than the profile in the original ObsSpace.
      const std::vector<std::size_t> &locsOriginal =
        data.obsspace().recidx_vector(recnums[jprof]);
      const std::vector<std::size_t> &locsExtended =
        data.obsspace().recidx_vector(recnums[jprof + nprofs]);

      // Indicates whether the averaged profile should be flagged.
      bool flagAveraged = false;
      // Flag the orginal profile in the same way as is done by the Reject action.
      for (size_t jloc : locsOriginal) {
        if (flagged[jv][jloc] && flags[iv][jloc] == QCflags::pass) {
          flags[iv][jloc] = filterQCflag;
          flagAveraged = true;
        }
      }

      // If any observation in the original profile was flagged with the filter QC flag,
      // flag all entries in the corresponding averaged profile.
      if (flagAveraged) {
        for (size_t jloc : locsExtended)
          flags[iv][jloc] = filterQCflag;
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
