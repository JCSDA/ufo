/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileMaxDifferenceCheck.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

ProfileMaxDifferenceCheck::ProfileMaxDifferenceCheck(
                                 ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ProfileMaxDifferenceCheck: "
                     << "using GNSSRO super refraction check method" << std::endl;
  allvars_ += parameters_.observationHeight;
  allvars_ += parameters_.variable;
  if (obsdb.obs_group_vars().empty())
    throw eckit::BadParameter("group variables configuration is empty.", Here());
}
// -----------------------------------------------------------------------------

ProfileMaxDifferenceCheck::~ProfileMaxDifferenceCheck() {
  oops::Log::trace() << "ProfileMaxDifferenceCheck: destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ProfileMaxDifferenceCheck::applyFilter(
                                      const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProfileMaxDifferenceCheck postFilter" << std::endl;
  const float missingFloat = util::missingValue<float>();

  // Get the variables that will be used in this filter
  std::vector<float> checkVariable;
  std::vector<float> observationHeight;
  data_.get(parameters_.variable, checkVariable);
  data_.get(parameters_.observationHeight, observationHeight);

  // Get the record numbers, so that we can group observations into profiles
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // Loop over the unique profiles
  for (size_t iProfile : record_numbers) {
    const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);
    // Find the set of indices which allow us to sort the variables by height
    std::vector<size_t> idx(obs_numbers.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
         [&observationHeight, &obs_numbers](size_t i1, size_t i2)
         {return observationHeight[obs_numbers[i1]] < observationHeight[obs_numbers[i2]];});

    // check if the profile does not penetrate below the maxCheckHeight
    if (observationHeight[obs_numbers[idx[0]]] >= parameters_.maxCheckHeight.value())
      continue;

    // Create the sorted heights and variables from the list of indices, checking
    // for missing values and whether we're within the check range
    std::vector<float> heightProfile;
    std::vector<float> variableProfile;
    for (size_t isort : idx) {
      if (observationHeight[obs_numbers[isort]] != missingFloat
          && checkVariable[obs_numbers[isort]] != missingFloat) {
        if (observationHeight[obs_numbers[isort]] <=
              (parameters_.maxCheckHeight.value()+parameters_.binSize.value())) {
          heightProfile.push_back(observationHeight[obs_numbers[isort]]);
          variableProfile.push_back(checkVariable[obs_numbers[isort]]);
        }
      }
    }
    //  Calculate the maximum different in the profile
    int jSuper;
    jSuper = calcMaxDifference(heightProfile, variableProfile);
    for (size_t iVar = 0; iVar < filtervars.nvars(); ++iVar) {
      for (size_t iobs : obs_numbers) {
        if (observationHeight[iobs] <= heightProfile[jSuper] &&
            jSuper != heightProfile.size()-1 )
          flagged[iVar][iobs] = true;
      }
    }  //  end iVar loop
  }  //  end iProfile loop
}  //  end applyFilter

int ProfileMaxDifferenceCheck::calcMaxDifference(
                                 const std::vector<float> & heightProfile,
                                 const std::vector<float> & variableProfile) const {
  int jSuper = heightProfile.size()-1;

  // Loop over all the heights within the range (starting at the top).
  // If the difference between the maximum and minimum values for a bin
  // above a certain height is above the threshold,
  // return the index where the threshold is first exceeded.
  for (size_t k = heightProfile.size()-1; k > 0; --k) {
    float minValue = 0;
    float maxValue = 0;
    if (heightProfile[k] > parameters_.maxCheckHeight.value()) {
      continue;
    } else {
      minValue = variableProfile[k];
      maxValue = variableProfile[k];
      for (size_t ibin = k+1; ibin < heightProfile.size(); ++ibin) {
        if (heightProfile[ibin] - heightProfile[k] >= parameters_.binSize.value()) {
          break;
        }
        minValue = std::min(minValue, variableProfile[ibin]);
        maxValue = std::max(maxValue, variableProfile[ibin]);
      }
      if ((maxValue - minValue) >= parameters_.threshold.value()) {
        jSuper = k;
        break;
      }
    }
  }  // end k loop
  return jSuper;
}
// -----------------------------------------------------------------------------

void ProfileMaxDifferenceCheck::print(std::ostream & os) const {
  os << "ProfileMaxDifferenceCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
