/*
 * (C) Crown copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/LocationToObsSpace.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"

namespace ufo {

static ObsFunctionMaker<LocationToObsSpace>
makerLocationToObsSpace_("LocationToObsSpace");

// -----------------------------------------------------------------------------

LocationToObsSpace::LocationToObsSpace(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "LocationToObsSpace constructor" << std::endl;
    options_.validateAndDeserialize(conf);
}

// -----------------------------------------------------------------------------

LocationToObsSpace::~LocationToObsSpace() {
  oops::Log::trace() << "LocationToObsSpace destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void LocationToObsSpace::compute(const ObsFilterData & in,
                                 ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "LocationToObsSpace compute start" << std::endl;

  const ioda::ObsSpace & obsdb = in.obsspace();

  if (options_.label_type.value() == LocationToObsSpaceLabelType::GlobalLocation) {
    const std::vector<std::size_t> & index = obsdb.index();
    for (std::size_t jloc = 0; jloc < obsdb.nlocs(); ++jloc) {
      out[0][jloc] = static_cast<int>(index[jloc]);
    }
  } else if (options_.label_type.value() == LocationToObsSpaceLabelType::InRecordLocation) {
    // Ensure observations have been grouped into records.
    if (obsdb.obs_group_vars().empty())
      throw eckit::UserError("Group variables configuration is empty", Here());

    const std::vector<size_t> & recordNumbers = obsdb.recidx_all_recnums();

    for (size_t iRecord : recordNumbers) {
      const std::vector<size_t> & recordIdxs = obsdb.recidx_vector(iRecord);
      for (size_t i = 0; i < recordIdxs.size(); ++i) {
        out[0][recordIdxs[i]] = static_cast<int>(i);
      }
    }
  }

  oops::Log::trace() << "LocationToObsSpace compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & LocationToObsSpace::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
