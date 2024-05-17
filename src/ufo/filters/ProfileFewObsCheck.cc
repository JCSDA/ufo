/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileFewObsCheck.h"

#include <algorithm>
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
  if ((!parameters_.threshold.value() && !parameters_.fraction.value()) ||
      (parameters_.threshold.value() && parameters_.fraction.value())) {
      throw eckit::BadParameter("Either threshold or fraction options must be specified.", Here());
  }
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
  const oops::ObsVariables observed = obsdb_.obsvariables();

  // Check the number of channels and variables to process
  const size_t nChans = std::max(obsdb_.nchans(), 1LU);
  // For multi-channel data "nvars" is the number of simulated variables times the
  // number of channels.  Therefore, divide "nvars" by the number of channels
  // to get the number of simulated variables.
  const size_t nActualVars = filtervars.nvars() / nChans;

  bool belowThreshold = false;

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // Loop over the number of actual variables
  for (size_t iVar=0; iVar < nActualVars; ++iVar) {
    std::vector<size_t> variableIndicesMap(nChans, -1);
    for (size_t iChan = 0; iChan < nChans; ++iChan) {
        const size_t iFilterVar = iVar * nChans + iChan;
        const std::string variableName = filtervars.variable(iFilterVar).variable();
        variableIndicesMap[iChan] = observed.find(variableName);
    }

    // Loop over the unique profiles
    for (size_t iProfile : record_numbers) {
      const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);

      int numValid = 0;
      int numTotal = 0;
      // For each channel and vertical level count the number of valid observations
      for (size_t iChan=0; iChan < nChans; ++iChan) {
        const size_t iFilterVar = iVar * nChans + iChan;
        const size_t jVar = variableIndicesMap[iChan];

        // Count the number of valid observations in this profile
        for (size_t jobs : obs_numbers) {
          if (apply[jobs]) {
            numTotal++;
            if ((*flags_)[jVar][jobs] == QCflags::pass)
              numValid++;
          }
        }
      }

      if (parameters_.threshold.value()) {
        belowThreshold = (numValid < parameters_.threshold.value().value());
      } else {
        belowThreshold = (static_cast<float>(numValid) / static_cast<float>(std::max(numTotal, 1)))
                         < parameters_.fraction.value().value();
      }

      oops::Log::debug() << "For var " << iVar << ", profile " << iProfile
                         << " there are " << numValid << "/" << numTotal << " valid observations"
                         << " belowThreshold = " << (belowThreshold ? "true" : "false")
                         << std::endl;
      if (belowThreshold) {
        // Reject profiles which don't contain sufficient observations
        for (size_t iChan=0; iChan < nChans; ++iChan) {
          // Note that this assumes that all the channels for each variable are
          // grouped together
          const size_t iFilterVar = iVar * nChans + iChan;
          const size_t jVar = variableIndicesMap[iChan];
          for (size_t jobs : obs_numbers)
            if (apply[jobs] && (*flags_)[jVar][jobs] == QCflags::pass)
              flagged[iFilterVar][jobs] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ProfileFewObsCheck::print(std::ostream & os) const {
  os << "ProfileFewObsCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
