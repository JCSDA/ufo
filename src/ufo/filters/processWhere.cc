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
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

ufo::Variables getAllWhereVariables(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  ufo::Variables vars;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    vars += ufo::Variables(masks[jm]);
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
  const int missing = util::missingValue(missing);
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (data[jj] == missing || oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

std::vector<bool> processWhere(const eckit::Configuration & config,
                               const ObsFilterData & filterdata) {
  const float missing = util::missingValue(missing);
  const size_t nlocs = filterdata.nlocs();

// Everywhere by default if no mask
  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  const ufo::Variables vars = getAllWhereVariables(config);
  for (size_t jm = 0; jm < vars.size(); ++jm) {
    if (vars.group(jm) != "VarMetaData") {
//    Process masks on float values
      const float vmin = masks[jm].getFloat("minvalue", missing);
      const float vmax = masks[jm].getFloat("maxvalue", missing);

//    Apply mask min/max
      if (vmin != missing || vmax != missing) {
        std::vector<float> data;
        filterdata.get(vars[jm], data);
        processWhereMinMax(data, vmin, vmax, where);
      }

//    Apply mask is_defined
      if (masks[jm].has("is_defined")) {
        if (filterdata.has(vars[jm])) {
          std::vector<float> data;
          filterdata.get(vars[jm], data);
          processWhereIsDefined(data, where);
        } else {
          std::fill(where.begin(), where.end(), false);
        }
      }

//    Apply mask is_not_defined
      if (masks[jm].has("is_not_defined")) {
        std::vector<float> data;
        filterdata.get(vars[jm], data);
        processWhereIsNotDefined(data, where);
      }

//    Apply mask is_in
      if (masks[jm].has("is_in")) {
        std::vector<int> data;
        filterdata.get(vars[jm], data);
        std::set<int> whitelist = oops::parseIntSet(masks[jm].getString("is_in"));
        processWhereIsIn(data, whitelist, where);
      }

//    Apply mask is_not_in
      if (masks[jm].has("is_not_in")) {
        std::vector<int> data;
        filterdata.get(vars[jm], data);
        std::set<int> blacklist = oops::parseIntSet(masks[jm].getString("is_not_in"));
        processWhereIsNotIn(data, blacklist, where);
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
