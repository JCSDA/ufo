/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/ReduceObsSpace.h"

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsSpace.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<ReduceObsSpace> makerReduceObsSpace_("reduce obs space");

// -----------------------------------------------------------------------------

ReduceObsSpace::ReduceObsSpace(const ReduceObsSpaceParameters &)
  : allvars_() {
}

// -----------------------------------------------------------------------------

void ReduceObsSpace::apply(const Variables &,
                           const std::vector<std::vector<bool>> & flagged,
                           ObsFilterData & filterdata,
                           int,
                           ioda::ObsDataVector<int> &,
                           ioda::ObsDataVector<float> &) const {
  if (filterdata.getGeoVaLs() != NULL) {
    throw eckit::NotImplemented("Action 'reduce obs space' not implemented for prior or post "
                                "filters yet");
  }
  const size_t nvars = flagged.size();
  if (nvars == 0) return;
  const size_t nlocs = flagged[0].size();
  std::vector<bool> keepObs(nlocs, false);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      // Keep observations if at least one of the variables needs to be kept
      if (!flagged[jvar][jloc]) keepObs[jloc] = true;
    }
  }
  filterdata.obsspace().reduce(keepObs);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
