/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/processWhere.h"

#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------

oops::Variables preProcessWhere(const eckit::Configuration & config, const std::string & group) {
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  oops::Variables vars;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    const std::string vargrp(masks[jm].getString("variable"));
    std::string var;
    std::string grp;
    splitVarGroup(vargrp, var, grp);
    if (grp == group) vars.push_back(var);
    if (group == "GeoVaLs" && grp.find("Function") != std::string::npos) {
      ObsFunction obsdiag(var);
      vars += obsdiag.requiredGeoVaLs();
    }
  }
  return vars;
}

// -----------------------------------------------------------------------------

void processWhereMinMax(const std::vector<float> & data,
            const float & vmin, const float & vmax, std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();

  if (vmin != missing || vmax != missing) {
    for (size_t jj = 0; jj < n; ++jj) {
      if (vmin != missing && data[jj] < vmin) mask[jj] = false;
      if (vmax != missing && data[jj] > vmax) mask[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------

void processWhereIsDefined(const std::vector<float> & data,
                           std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (data[jj] == missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsNotDefined(const std::vector<float> & data,
                              std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (data[jj] != missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsIn(const std::vector<int> & data,
                      const std::set<int> & whitelist,
                      std::vector<bool> & mask) {
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (!oops::contains(whitelist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsNotIn(const std::vector<int> & data,
                         const std::set<int> & blacklist,
                         std::vector<bool> & mask) {
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

std::vector<bool> processWhere(const eckit::Configuration & config,
                               ObsFilterData & filterdata) {
  const float missing = util::missingValue(missing);
  const size_t nlocs = filterdata.nlocs();

// Everywhere by default if no mask
  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  for (size_t jm = 0; jm < masks.size(); ++jm) {
//  Get variable@group
    const std::string varname(masks[jm].getString("variable"));
    std::string var, grp;
    splitVarGroup(varname, var, grp);
    if (grp != "VarMetaData") {
//    Get data
      std::vector<float> data = filterdata.get(varname);
//    Process masks on float values
      const float vmin = masks[jm].getFloat("minvalue", missing);
      const float vmax = masks[jm].getFloat("maxvalue", missing);

//    Apply mask min/max
      if (vmin != missing || vmax != missing) {
        processWhereMinMax(data, vmin, vmax, where);
      }

//    Apply mask is_defined
      if (masks[jm].has("is_defined")) {
        if (filterdata.has(varname)) {
          processWhereIsDefined(data, where);
        } else {
          std::fill(where.begin(), where.end(), false);
        }
      }

//    Apply mask is_not_defined
      if (masks[jm].has("is_not_defined")) {
        processWhereIsNotDefined(data, where);
      }
    }
  }
//  Print diagnostics for debug
  int ii = 0;
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (where[jj] == false) ++ii;
  }

  oops::Log::debug() << "processWhere: selected " << ii << " obs." << std::endl;
  return where;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
