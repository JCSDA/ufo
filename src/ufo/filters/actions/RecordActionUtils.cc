/*
 * (C) Copyright 2026 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/RecordActionUtils.h"

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsSpace.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {
namespace actions {

std::vector<std::vector<std::size_t>> recordLocationsOrThrow(
    const ObsFilterData &data, const std::string &actionName) {
  // Check that the ObsSpace is divided into records (groups).
  // If not, record-mode actions cannot proceed.
  if (data.obsspace().obs_group_vars().empty()) {
    throw eckit::BadParameter(
        "The ObsSpace must have been divided into records to use the '" +
            actionName + "' action with 'apply to whole record: true', but "
            "the group variables configuration is empty.",
        Here());
  }

  // Get all record numbers from the ObsSpace.
  const std::vector<std::size_t> recnums = data.obsspace().recidx_all_recnums();

  // For each record number, retrieve the vector of observation indices in that record.
  std::vector<std::vector<std::size_t>> recordLocs;
  recordLocs.reserve(recnums.size());
  for (const std::size_t recnum : recnums) {
    recordLocs.push_back(data.obsspace().recidx_vector(recnum));
  }

  return recordLocs;
}

}  // namespace actions
}  // namespace ufo
