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

#include "ufo/filters/QCflags.h"
#include "ufo/utils/RecordHandler.h"

namespace ufo {
  RecordHandler::RecordHandler(const ioda::ObsSpace & obsdb,
                               const Variables & filtervars,
                               const ioda::ObsDataVector<int> & flags,
                               const bool retainOnlyIfAllFilterVariablesAreValid)
    : obsdb_(obsdb),
      filtervars_(filtervars),
      flags_(flags),
      retainOnlyIfAllFilterVariablesAreValid_(retainOnlyIfAllFilterVariablesAreValid)
  {}

std::vector<std::size_t> RecordHandler::getLaunchPositions() const {
  const util::DateTime missingDateTime = util::missingValue<util::DateTime>();

  // Retrieve datetimes.
  std::vector<util::DateTime> dateTimes(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "dateTime", dateTimes);

  // Indices of filter variables in the `flags` vector.
  std::vector<size_t> indexOfFilterVariableInFlags;
  for (size_t ivar = 0; ivar < filtervars_.nvars(); ++ivar) {
    const std::string filterVariableName = filtervars_.variable(ivar).variable();
    indexOfFilterVariableInFlags.push_back(flags_.varnames().find(filterVariableName));
  }

  // Vector of locations corresponding to profile launch positions.
  std::vector<std::size_t> launchPositions;

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = obsdb_.recidx_all_recnums();
  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < recnums.size(); ++jprof) {
    // Get locations corresponding to this profile.
    std::vector<std::size_t> locs = obsdb_.recidx_vector(recnums[jprof]);
    // Sort locs according to values of dateTime, ignoring missing values.
    std::stable_sort(locs.begin(),
                     locs.end(),
                     [&](int a, int b)
                     {return dateTimes[a] == missingDateTime ?
                         false :
                         (dateTimes[b] == missingDateTime ?
                          true :
                          dateTimes[a] < dateTimes[b]);});

    // Find the location corresponding to the launch position.
    // This is defined as the location with the earliest non-missing datetime
    // and a certain number of filter variables with QC flags equal to pass.
    // If `retainOnlyIfAllFilterVariablesAreValid` is true, all filter variables must
    // have QC flags equal to pass. If it is false then at least one must have a
    // QC flag equal to pass.
    size_t launchPosition = locs.front();
    for (const size_t jloc : locs) {
      // Skip location if dateTime is missing.
      if (dateTimes[jloc] == missingDateTime) continue;
      // Skip location if a certain number of filter variables are missing.
      bool filterVarsOK = true;
      if (retainOnlyIfAllFilterVariablesAreValid_) {
        for (const int idx : indexOfFilterVariableInFlags) {
          if (QCflags::isRejected(flags_[idx][jloc])) {
            filterVarsOK = false;
            break;
          }
        }
      } else {
        filterVarsOK = false;
        for (const int idx : indexOfFilterVariableInFlags) {
          if (!QCflags::isRejected(flags_[idx][jloc])) {
            filterVarsOK = true;
            break;
          }
        }
      }
      if (!filterVarsOK) continue;
      // Launch position is current location.
      launchPosition = jloc;
      break;
    }
    launchPositions.push_back(launchPosition);
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

