/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FilterUtils.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

namespace {

bool shouldUnselectLocation(size_t location,
                            const std::vector<ioda::ObsDataRow<int>> &flags,
                            UnselectLocationIf mode) {
  switch (mode) {
  case UnselectLocationIf::ANY_FILTER_VARIABLE_REJECTED:
    for (size_t variable = 0; variable < flags.size(); ++variable)
      if (QCflags::isRejected(flags[variable][location]))
        return true;
    return false;

  case UnselectLocationIf::ALL_FILTER_VARIABLES_REJECTED:
    for (size_t variable = 0; variable < flags.size(); ++variable)
      if (!QCflags::isRejected(flags[variable][location]))
        return false;
    return true;
  }

  // Should not get here
  return false;
}

}  // namespace

void unselectRejectedLocations(std::vector<bool> &selected,
                               const ufo::Variables &filtervars,
                               const ioda::ObsDataVector<int> &qcflags,
                               UnselectLocationIf mode,
                               const std::vector<size_t> &obs_inds) {
  std::vector<ioda::ObsDataRow<int>> filterVariableFlags;
  // Select flags for respective filtervars
  for (size_t ivar = 0; ivar < filtervars.nvars(); ++ivar) {
    const std::string filterVariableName = filtervars.variable(ivar).variable();
    filterVariableFlags.push_back(qcflags[filterVariableName]);
  }
  const bool specificInds = (obs_inds.size() > 0);
  const size_t nLoc = specificInds ? obs_inds.size() : selected.size();
  for (size_t iloc = 0; iloc < nLoc; ++iloc) {
    const size_t loc = specificInds ? obs_inds[iloc] : iloc;
    if (selected[loc] && shouldUnselectLocation(loc, filterVariableFlags, mode))
      selected[loc] = false;
  }
}


}  // namespace ufo
