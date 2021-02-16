/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileFewObsCheck.h"

#include <iomanip>
#include <iostream>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// ProfileFewObsCheck: Check the number of observations in a profile
///
/// This will use the record number in obsdb_ to identify which observations belong to a
/// given profile (all members of a profile must share the same record number).
/// For each profile the number of valid observations is found, and the profile
/// is rejected if this is below the given threshold.

ProfileFewObsCheck::ProfileFewObsCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ProfileFewObsCheck constructor" << std::endl;
}

// -----------------------------------------------------------------------------

ProfileFewObsCheck::~ProfileFewObsCheck() {
  oops::Log::trace() << "ProfileFewObsCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the profile check for the number of observations.

void ProfileFewObsCheck::applyFilter(const std::vector<bool> & apply,
                                     const Variables & filtervars,
                                     std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProfileFewObsCheck preProcess filter" << std::endl;
  const oops::Variables observed = obsdb_.obsvariables();

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // For each variable, check the number of observations in the profile
  for (size_t iFilterVar = 0; iFilterVar < filtervars.nvars(); ++iFilterVar) {
    const size_t iVar = observed.find(filtervars.variable(iFilterVar).variable());

    // Loop over the unique profiles
    for (size_t iProfile : record_numbers) {
      const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);

      // Count the number of valid observations in this profile
      int numValid = 0;
      for (size_t jobs : obs_numbers)
        if (apply[jobs] && (*flags_)[iVar][jobs] == QCflags::pass)
          numValid++;

      // Reject profiles which don't contain sufficient observations
      if (numValid < parameters_.threshold.value())
        for (size_t jobs : obs_numbers)
          if (apply[jobs] && (*flags_)[iVar][jobs] == QCflags::pass)
            flagged[iFilterVar][jobs] = true;
    }
  }
}

// -----------------------------------------------------------------------------

void ProfileFewObsCheck::print(std::ostream & os) const {
  os << "<FilterName>: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
