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
#include "oops/base/Variables.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

oops::Variables preProcessWhere(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  std::vector<std::string> vv;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    const std::string vargrp(masks[jm].getString("variable"));
    std::string var;
    std::string grp;
    splitVarGroup(vargrp, var, grp);
    if (grp == "GeoVaLs") vv.push_back(var);
  }

  return oops::Variables(vv);
}

// -----------------------------------------------------------------------------

std::vector<bool> processWhere(ioda::ObsSpace & obsdb, const GeoVaLs & gvals,
                               const eckit::Configuration & config) {
  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb.nlocs();

// Everywhere by default if no mask
  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  for (size_t jm = 0; jm < masks.size(); ++jm) {
//  Get variable and group
    const std::string vargrp(masks[jm].getString("variable"));
    std::string var = vargrp;
    std::string grp = "MetaData";
    splitVarGroup(vargrp, var, grp);

//  Set obs group if group is not GeoVaLs
    std::string obgrp = grp;
    if (grp == "GeoVaLs") obgrp = "";
    if (obgrp == "MetaData") {                               // TEMPORARY HACK TO BE REMOVED
      if (!obsdb.has(obgrp, var)) obgrp = "GroupUndefined";  // TEMPORARY HACK TO BE REMOVED
    }                                                        // TEMPORARY HACK TO BE REMOVED

//  Process masks on float values
    const float vmin = masks[jm].getFloat("minvalue", missing);
    const float vmax = masks[jm].getFloat("maxvalue", missing);

    if (vmin != missing || vmax != missing ||
        masks[jm].has("is_defined") || masks[jm].has("is_not_defined")) {
//    Get float values
      ioda::ObsDataVector<float> values(obsdb, var, obgrp);
      if (grp == "GeoVaLs") gvals.get(values.values(), var);
      if (grp == "GeoVaLs") oops::Log::debug() << "processWhere gv = " << gvals << std::endl;
      if (grp == "GeoVaLs") oops::Log::debug() << "processWhere values = " << values << std::endl;

//    Apply mask min/max
      if (vmin != missing || vmax != missing) {
        for (size_t jj = 0; jj < nlocs; ++jj) {
          if (vmin != missing && values[jj] < vmin) where[jj] = false;
          if (vmax != missing && values[jj] > vmax) where[jj] = false;
        }
      }

//    Apply mask is_defined
      if (masks[jm].has("is_defined")) {
        if (obsdb.has(obgrp, var) || grp == "GeoVaLs") {
          for (size_t jj = 0; jj < nlocs; ++jj) {
            if (values[jj] == missing) where[jj] = false;
          }
        } else {
          for (size_t jj = 0; jj < nlocs; ++jj) where[jj] = false;
        }
      }

//    Apply mask is_not_defined
      if (masks[jm].has("is_not_defined")) {
        if (obsdb.has(obgrp, var) || grp == "GeoVaLs") {
          for (size_t jj = 0; jj < nlocs; ++jj) {
            if (values[jj] != missing) where[jj] = false;
          }
        }
      }
    }

//  Process masks on integer values
    if (masks[jm].has("is_in") || masks[jm].has("is_not_in")) {
//    Get int values
      ioda::ObsDataVector<int> valint(obsdb, var, obgrp);

//    Apply mask is_in
      if (masks[jm].has("is_in")) {
        std::set<int> whitelist = parseIntSet(masks[jm].getString("is_in"));
        for (size_t jj = 0; jj < nlocs; ++jj) {
          if (!contains(whitelist, valint[jj])) where[jj] = false;
        }
      }

//    Apply mask is_not_in
      if (masks[jm].has("is_not_in")) {
        std::set<int> blacklist = parseIntSet(masks[jm].getString("is_not_in"));
        for (size_t jj = 0; jj < nlocs; ++jj) {
          if (contains(blacklist, valint[jj])) where[jj] = false;
        }
      }
    }
  }

// Print diagnostics for debug
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
