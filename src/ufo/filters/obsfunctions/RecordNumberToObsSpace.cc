/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/RecordNumberToObsSpace.h"

#include <vector>

#include "ioda/ObsDataVector.h"

namespace ufo {

static ObsFunctionMaker<RecordNumberToObsSpace>
makerRecordNumberToObsSpace_("RecordNumberToObsSpace");

// -----------------------------------------------------------------------------

RecordNumberToObsSpace::RecordNumberToObsSpace(const eckit::LocalConfiguration & conf)
  : invars_() {
}

// -----------------------------------------------------------------------------

RecordNumberToObsSpace::~RecordNumberToObsSpace() {}

// -----------------------------------------------------------------------------

void RecordNumberToObsSpace::compute(const ObsFilterData & in,
                                     ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "RecordNumberToObsSpace::compute started" << std::endl;

  const ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into records.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Correspondence between record numbers and locations in the data sample.
  const std::vector<std::size_t> & recnums = obsdb.recnum();

  for (std::size_t jloc = 0; jloc < obsdb.nlocs(); ++jloc) {
    out[0][jloc] = static_cast<int>(recnums[jloc]);
  }

  oops::Log::trace() << "RecordNumberToObsSpace::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & RecordNumberToObsSpace::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
