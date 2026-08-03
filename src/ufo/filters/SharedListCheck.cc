/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SharedListCheck.h"

#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"

#include "ufo/filters/SharedListStore.h"
#include "ufo/filters/Variable.h"

namespace ufo {

// -----------------------------------------------------------------------------

SharedListCheck::SharedListCheck(ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 ioda::ObsDataVector<int> & flags,
                                 ioda::ObsDataVector<float> & obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::trace() << "SharedListCheck constructor" << std::endl;

  // Load the shared list file (lazy — skips if already loaded)
  SharedListStore::instance().load(parameters_.sharedListFile.value());
}

// -----------------------------------------------------------------------------

SharedListCheck::~SharedListCheck() {
  oops::Log::trace() << "SharedListCheck destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void SharedListCheck::applyFilter(
    const std::vector<bool> & apply,
    const Variables & filtervars,
    std::vector<std::vector<bool>> & flagged) const {

  oops::Log::trace() << "SharedListCheck applyFilter start" << std::endl;

  const std::string & filepath = parameters_.sharedListFile.value();
  const std::string & mode = parameters_.flagMode.value();
  const eckit::LocalConfiguration & conf =
      SharedListStore::instance().getConfig(filepath);

  const std::string & listKey = parameters_.listToUse.value();

  // Determine matching mode: compound check vs simple
  const bool isCompound =
      (parameters_.compoundCheck.value() != boost::none);

  // Loop over filter variables
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    // Skip if this key is not in the list file
    if (!conf.has(listKey)) {
      continue;
    }

    if (isCompound) {
      // Compound matching: match multiple obs fields against corresponding lists at the same time
      const auto & checkEntries = *parameters_.compoundCheck.value();

      // Read all obs variables for compound check
      std::vector<std::vector<std::string>> obsValues(checkEntries.size());
      std::vector<std::string> sublistsToUse(checkEntries.size());
      for (size_t ic = 0; ic < checkEntries.size(); ++ic) {
        const Variable obsVar(checkEntries[ic].obsVariable.value());
        obsValues[ic].resize(obsdb_.nlocs());
        obsdb_.get_db(obsVar.group(), obsVar.variable(), obsValues[ic]);
        sublistsToUse[ic] = checkEntries[ic].sublistToUse.value();
      }

      // Read list entries as vector of configurations
      std::vector<eckit::LocalConfiguration> entries;
      conf.get(listKey, entries);

      for (size_t loc = 0; loc < obsdb_.nlocs(); ++loc) {
        if (!apply[loc]) continue;

        bool matched = false;
        for (const auto & entry : entries) {
          bool entryMatch = true;
          for (size_t ic = 0; ic < checkEntries.size(); ++ic) {
            std::vector<std::string> listVals;
            entry.get(sublistsToUse[ic], listVals);
            // Check if "all" is in the list or obs value matches any entry
            bool fieldMatch = false;
            for (const auto & val : listVals) {
              if (val == "all" || val == obsValues[ic][loc]) {
                fieldMatch = true;
                break;
              }
            }
            if (!fieldMatch) {
              entryMatch = false;
              break;
            }
          }
          if (entryMatch) {
            matched = true;
            break;
          }
        }

        if (mode == "flag matched") {
          flagged[jv][loc] = matched;
        } else if (mode == "flag unmatched") {
          flagged[jv][loc] = !matched;
        }
      }
    } else {
      // Simple matching (single variable lookup in set)
      const Variable obsKeyVar(*parameters_.variableToCheck.value());
      std::vector<std::string> obsKeys(obsdb_.nlocs());
      obsdb_.get_db(obsKeyVar.group(), obsKeyVar.variable(), obsKeys);

      std::vector<std::string> ids;
      conf.get(listKey, ids);
      const std::set<std::string> idSet(ids.begin(), ids.end());

      for (size_t loc = 0; loc < obsdb_.nlocs(); ++loc) {
        if (!apply[loc]) continue;

        const bool matched = (idSet.count(obsKeys[loc]) > 0);

        if (mode == "flag matched") {
          flagged[jv][loc] = matched;
        } else if (mode == "flag unmatched") {
          flagged[jv][loc] = !matched;
        }
      }
    }
  }

  oops::Log::trace() << "SharedListCheck applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void SharedListCheck::print(std::ostream & os) const {
  os << "SharedListCheck: file=" << parameters_.sharedListFile.value()
     << ", flag mode=" << parameters_.flagMode.value();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
