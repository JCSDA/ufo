/*
 * (C) Copyright 2026 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_RECORDACTIONUTILS_H_
#define UFO_FILTERS_ACTIONS_RECORDACTIONUTILS_H_

#include <cstddef>
#include <string>
#include <vector>

namespace ufo {

class ObsFilterData;

namespace actions {

/// \brief Retrieve the observation locations for each record from ObsSpace.
///
/// Throws eckit::BadParameter if the ObsSpace has not been divided into records
/// (i.e., if obs_group_vars() is empty).
/// \param data Accessor to obs filter data containing the ObsSpace.
/// \param actionName Name of the filter action (used in error message).
/// \return Vector of vectors; recordLocs[i] contains the observation indices
/// for record i.
std::vector<std::vector<std::size_t>> recordLocationsOrThrow(
    const ObsFilterData &data, const std::string &actionName);

/// \brief Expand per-location flagged mask to whole-record flagged mask.
///
/// This is the core of record-mode action processing. It:
///   1. Creates an expanded mask initialized to false.
///   2. For each variable and each record, checks if any location in that
///   record
///      is (a) flagged and (b) satisfies the provided eligibility check.
///   3. If yes, marks the entire record as true for that variable.
///
/// Each action supplies a custom eligibility check that encodes the specific
/// per-location eligibility criteria used in its apply() method.
///
/// Template parameter Predicate:
///   A callable that accepts (ifiltervar, jloc) and returns bool.
///   Examples:
///     - AcceptObs: check that the current QC flag can be accepted.
///     - RejectObs: also check if the current QC is "pass".
///     - SetFlag: also check if the observation is not "ignored".
///
/// \param flagged 2D bool array of (variable, location) flags. If true, this
///        location should trigger a whole record update as long as the
///        eligibility check is satisfied.
/// \param recordLocs Vector of vectors; recordLocs[i] = vector of location
/// indices in record i.
/// \param isEligibleForRecordExpansion Callable(ifiltervar, jloc) -> bool;
/// returns true if this
///        (variable, location)
///        pair is eligible to trigger record expansion.
/// \return A new 2D bool array (same shape as \p flagged) where each variable's
/// mask
///         is expanded so that entire records are marked true if any location
///         in that record was flagged and satisfied the eligibility check.
template <typename Predicate>
std::vector<std::vector<bool>> expandFlaggedToWholeRecord(
    const std::vector<std::vector<bool>> &flagged,
    const std::vector<std::vector<std::size_t>> &recordLocs,
    const Predicate &isEligibleForRecordExpansion) {
  // Initialize expanded mask to same shape as input, all false.
  std::vector<std::vector<bool>> expandedFlagged;
  expandedFlagged.reserve(flagged.size());
  for (const auto &flaggedVar : flagged) {
    expandedFlagged.emplace_back(flaggedVar.size(), false);
  }

  // For each filter variable and each record:
  for (std::size_t ifiltervar = 0; ifiltervar < flagged.size(); ++ifiltervar) {
    for (const auto &locsInRecord : recordLocs) {
      // Check if any location in this record is flagged and satisfies the
      // eligibility check.
      bool anyFlaggedInThisRecord = false;
      for (std::size_t jloc : locsInRecord) {
        if (flagged[ifiltervar][jloc] &&
            isEligibleForRecordExpansion(ifiltervar, jloc)) {
          anyFlaggedInThisRecord = true;
          break;
        }
      }

      if (anyFlaggedInThisRecord) {
        // Mark all locations in this record as true for this variable.
        for (std::size_t jloc : locsInRecord) {
          expandedFlagged[ifiltervar][jloc] = true;
        }
      }
    }
  }

  return expandedFlagged;
}

}  // namespace actions
}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_RECORDACTIONUTILS_H_
