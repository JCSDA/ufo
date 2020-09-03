/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/processWhere.h"

#include <bitset>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

ufo::Variables getAllWhereVariables(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  ufo::Variables vars;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    eckit::LocalConfiguration varconf(masks[jm], "variable");
    vars += ufo::Variable(varconf);
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
template <class T>
void processWhereIsIn(const std::vector<T> & data,
                      const std::set<T> & whitelist,
                      std::vector<bool> & mask) {
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (!oops::contains(whitelist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
template <class T>
void processWhereIsNotIn(const std::vector<T> & data,
                         const std::set<T> & blacklist,
                         std::vector<bool> & mask) {
  const T missing = util::missingValue(missing);
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (data[jj] == missing || oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void processWhereIsNotIn(const std::vector<std::string> & data,
                         const std::set<std::string> & blacklist,
                         std::vector<bool> & mask) {
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void processWhereBitSet(const std::vector<int> & data,
                        const std::set<int> & flags,
                        std::vector<bool> & mask) {
  std::bitset<32> flags_bs;
  for (const int &elem : flags) {
    flags_bs[elem] = 1;
  }
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if ((data[jj] & flags_bs.to_ulong()) != 0) mask[jj] = false;
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

  for (size_t jm = 0; jm < masks.size(); ++jm) {
    eckit::LocalConfiguration varconf(masks[jm], "variable");
    Variable var(varconf);
    for (size_t jvar = 0; jvar < var.size(); ++jvar) {
      if (var.group() != "VarMetaData") {
        const Variable varname = var[jvar];
        ioda::ObsDtype dtype = filterdata.dtype(varname);

//      Process masks on float values
        const float vmin = masks[jm].getFloat("minvalue", missing);
        const float vmax = masks[jm].getFloat("maxvalue", missing);

//      Apply mask min/max
        if (vmin != missing || vmax != missing) {
          std::vector<float> data;
          filterdata.get(varname, data);
          processWhereMinMax(data, vmin, vmax, where);
        }

//      Apply mask is_defined
        if (masks[jm].has("is_defined")) {
          if (filterdata.has(varname)) {
            std::vector<float> data;
            filterdata.get(varname, data);
            processWhereIsDefined(data, where);
          } else {
            std::fill(where.begin(), where.end(), false);
          }
        }

//      Apply mask is_not_defined
        if (masks[jm].has("is_not_defined")) {
          std::vector<float> data;
          filterdata.get(varname, data);
          processWhereIsNotDefined(data, where);
        }

//      Apply mask is_in
        if (masks[jm].has("is_in")) {
          if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            std::vector<std::string> whitelistvec = masks[jm].getStringVector("is_in");
            std::set<std::string> whitelist(whitelistvec.begin(), whitelistvec.end());
            filterdata.get(varname, data);
            processWhereIsIn(data, whitelist, where);
          } else if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            std::set<int> whitelist = oops::parseIntSet(masks[jm].getString("is_in"));
            filterdata.get(varname, data);
            processWhereIsIn(data, whitelist, where);
          } else {
            throw eckit::UserError(
              "Only integer and string variables may be used for processWhere 'is_in'",
              Here());
          }
        }

//      Apply mask is_not_in
        if (masks[jm].has("is_not_in")) {
          if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            std::vector<std::string> blacklistvec = masks[jm].getStringVector("is_not_in");
            std::set<std::string> blacklist(blacklistvec.begin(), blacklistvec.end());
            filterdata.get(varname, data);
            processWhereIsNotIn(data, blacklist, where);
          } else if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            std::set<int> blacklist = oops::parseIntSet(masks[jm].getString("is_not_in"));
            processWhereIsNotIn(data, blacklist, where);
          } else {
            throw eckit::UserError(
              "Only integer and string variables may be used for processWhere 'is_not_in'",
              Here());
          }
        }

//      Apply mask any_bit_set_of
        if (masks[jm].has("any_bit_set_of")) {
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            std::set<int> flags = oops::parseIntSet(masks[jm].getString("any_bit_set_of"));
            filterdata.get(varname, data);
            processWhereBitSet(data, flags, where);
          } else {
            throw eckit::UserError(
              "Only integer variables may be used for processWhere 'any_bit_set_of'",
              Here());
          }
        }
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
