/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/processWhere.h"

#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/missingValues.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

std::vector<bool> processWhere(ioda::ObsSpace & obsdb, const eckit::Configuration & config) {
  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb.nlocs();

  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  for (size_t jm = 0; jm < masks.size(); ++jm) {
    const std::string var(masks[jm].getString("variable"));
    std::string obgrp = "MetaData";
    if (!obsdb.has(obgrp, var)) obgrp = "GroupUndefined";  // TO BE REMOVED

    const float vmin = masks[jm].getFloat("minvalue", missing);
    const float vmax = masks[jm].getFloat("maxvalue", missing);

    if (vmin != missing || vmax != missing) {
      ioda::ObsDataVector<float> values(obsdb, var, obgrp);
      for (size_t jj = 0; jj < nlocs; ++jj) {
        if (vmin != missing && values[jj] < vmin) where[jj] = false;
        if (vmax != missing && values[jj] > vmax) where[jj] = false;
      }
    }

    if (masks[jm].has("is_in")) {
      ioda::ObsDataVector<int> values(obsdb, var, obgrp);
      std::set<int> whitelist = parseIntSet(masks[jm].getString("is_in"));
      for (size_t jj = 0; jj < nlocs; ++jj) {
        if (!contains(whitelist, values[jj])) where[jj] = false;
      }
    }

    if (masks[jm].has("is_not_in")) {
      ioda::ObsDataVector<int> values(obsdb, var, obgrp);
      std::set<int> blacklist = parseIntSet(masks[jm].getString("is_not_in"));
      for (size_t jj = 0; jj < nlocs; ++jj) {
        if (contains(blacklist, values[jj])) where[jj] = false;
      }
    }
  }

  int ii = 0;
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (where[jj] == false) ++ii;
  }
  oops::Log::debug() << "processWhere: " << obsdb.obsname()
                     << " selected " << ii << " obs." << std::endl;
  return where;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
