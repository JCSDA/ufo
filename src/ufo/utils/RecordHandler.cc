/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>

#include "ioda/ObsSpace.h"

#include "oops/util/missingValues.h"

#include "ufo/utils/RecordHandler.h"

namespace ufo {
  RecordHandler::RecordHandler(const ioda::ObsSpace & obsdb)
    : obsdb_(obsdb)
  {}

std::vector<std::size_t> RecordHandler::getLaunchPositions() const {
  const util::DateTime missingDateTime = util::missingValue(missingDateTime);

  // Retrieve datetimes.
  std::vector<util::DateTime> dateTimes(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "dateTime", dateTimes);

  // Vector of locations corresponding to profile launch positions.
  std::vector<std::size_t> launchPositions;

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = obsdb_.recidx_all_recnums();
  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < recnums.size(); ++jprof) {
    // Get locations corresponding to this profile.
    const std::vector<std::size_t> & locs = obsdb_.recidx_vector(recnums[jprof]);
    // Find the location corresponding to the launch position.
    // This is defined as the smallest non-missing datetime in the profile.
    // If all datetimes are missing this will select the first entry in the profile.
    auto it_nonmissing = std::min_element(locs.begin(), locs.end(),
                                          [missingDateTime, & dateTimes]
                                          (const size_t & a, const size_t & b)
                                          {return dateTimes[a] == missingDateTime ?
                                           false :
                                           (dateTimes[b] == missingDateTime ?
                                            true :
                                            dateTimes[a] < dateTimes[b]);});
    launchPositions.push_back(*it_nonmissing);
  }

  return launchPositions;
}

std::vector<bool> RecordHandler::changeApplyIfRecordsAreSingleObs
(const std::vector<bool> & apply) const {
  // The ObsSpace must have been grouped into records.
  // If not, the input vector is returned without modification.
  if (obsdb_.obs_group_vars().empty())
    return apply;

  // Create a vector for which the value of apply at the launch position
  // of each record is set to the logical `or` of the apply values in
  // that record. All other values in the vector are set to false.
  std::vector<bool> applyRecord(obsdb_.nlocs(), false);

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = obsdb_.recidx_all_recnums();

  // Loop over profiles.
  const std::vector<std::size_t> launchPositions = getLaunchPositions();
  for (std::size_t jprof = 0; jprof < recnums.size(); ++jprof) {
    // Get locations corresponding to this profile.
    const std::vector<std::size_t> & locs = obsdb_.recidx_vector(recnums[jprof]);
    // Logical `or` of values of apply in the entire record.
    bool applyLogicalOr = false;
    for (std::size_t loc : locs) {
      applyLogicalOr = applyLogicalOr || apply[loc];
      if (applyLogicalOr) break;
    }
    // Set apply at launch position to logical `or` of apply values.
    applyRecord[launchPositions[jprof]] = applyLogicalOr;
  }

  return applyRecord;
}

std::vector<bool> RecordHandler::changeThinnedIfRecordsAreSingleObs
(const std::vector<bool> & isThinned) const {
  // The ObsSpace must have been grouped into records.
  // If not, the input vector is returned without modification.
  if (obsdb_.obs_group_vars().empty())
    return isThinned;

  // Create a vector for which the values of isThinned at all locations
  // in each record are set to logical `or` of the isThinned values in that record.
  std::vector<bool> isThinnedRecord(isThinned.size(), false);

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = obsdb_.recidx_all_recnums();

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < recnums.size(); ++jprof) {
    // Get locations corresponding to this profile.
    const std::vector<std::size_t> & locs = obsdb_.recidx_vector(recnums[jprof]);
    // Logical `or` of values of isThinned in the entire record.
    bool isThinnedLogicalOr = false;
    for (std::size_t loc : locs) {
      const std::size_t gloc = obsdb_.distribution()->globalUniqueConsecutiveLocationIndex(loc);
      isThinnedLogicalOr = isThinnedLogicalOr || isThinned[gloc];
      if (isThinnedLogicalOr) break;
    }
    // Set isThinned at all locations in record to logical `or` of isThinned values.
    for (std::size_t loc : locs) {
      const std::size_t gloc = obsdb_.distribution()->globalUniqueConsecutiveLocationIndex(loc);
      isThinnedRecord[gloc] = isThinnedLogicalOr;
    }
  }

  return isThinnedRecord;
}

void RecordHandler::checkRecordCategories(const Variable & categoryVariableName) const {
  switch (obsdb_.dtype(categoryVariableName.group(),
                       categoryVariableName.variable())) {
  case ioda::ObsDtype::Integer:
    checkRecordCategoriesImpl<int>(categoryVariableName);
    break;
  case ioda::ObsDtype::String:
    checkRecordCategoriesImpl<std::string>(categoryVariableName);
    break;
  default:
    throw eckit::UserError(categoryVariableName.fullName() +
          " is neither an integer nor a string variable", Here());
  }
}

template <typename VariableType>
void RecordHandler::checkRecordCategoriesImpl(const Variable & categoryVariableName) const {
  // Get category variable from ObsSpace.
  std::vector <VariableType> categoryVariable(obsdb_.nlocs());
  obsdb_.get_db(categoryVariableName.group(),
                categoryVariableName.variable(),
                categoryVariable);

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = obsdb_.recidx_all_recnums();

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < recnums.size(); ++jprof) {
    // Get locations corresponding to this profile.
    const std::vector<std::size_t> & locs = obsdb_.recidx_vector(recnums[jprof]);
    for (std::size_t loc : locs) {
      if (categoryVariable[loc] != categoryVariable[locs.front()]) {
        throw eckit::UserError("Cannot have multiple categories per record", Here());
      }
    }
  }
}

}  // namespace ufo

