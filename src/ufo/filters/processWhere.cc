/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/processWhere.h"

#include <bitset>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PartialDateTime.h"
#include "oops/util/wildcard.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

namespace ufo {


// -----------------------------------------------------------------------------
ufo::Variables getAllWhereVariables(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> masks = config.getSubConfigurations();

  ufo::Variables vars;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    eckit::LocalConfiguration varconf(masks[jm], "variable");
    vars += ufo::Variable(varconf);
  }
  return vars;
}


// -----------------------------------------------------------------------------
template<typename T>
void processWhereMinMax(const std::vector<T> & data,
                        const T & vmin, const T & vmax,
                        std::vector<bool> & mask) {
  const T not_set_value = util::missingValue(not_set_value);
  const size_t n = data.size();

  if (vmin != not_set_value || vmax != not_set_value) {
    for (size_t jj = 0; jj < n; ++jj) {
      if (vmin != not_set_value && data[jj] < vmin) mask[jj] = false;
      if (vmax != not_set_value && data[jj] > vmax) mask[jj] = false;
    }
  }
}


// -----------------------------------------------------------------------------
void processWhereMinMax(const std::vector<util::DateTime> & data,
                        const std::string & vmin, const std::string & vmax,
                        std::vector<bool> & mask) {
  const std::string not_set_value = "0000-00-00T00:00:00Z";

  if (vmin != not_set_value || vmax != not_set_value) {
    util::PartialDateTime pdt_vmin(vmin), pdt_vmax(vmax);

    for (size_t jj = 0; jj < data.size(); ++jj) {
      if (vmin != not_set_value && pdt_vmin > data[jj]) mask[jj] = false;
      if (vmax != not_set_value && pdt_vmax < data[jj]) mask[jj] = false;
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
template <typename T>
void applyMinMax(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                 ObsFilterData const & filterdata, Variable const & varname) {
  const T not_set_value = util::missingValue(not_set_value);

  T vmin = not_set_value;
  // Set vmin to the value of the 'minvalue' option if it exists; if not, leave vmin unchanged.
  mask.get("minvalue", vmin);
  T vmax = not_set_value;
  mask.get("maxvalue", vmax);

  // Apply mask min/max
  if (vmin != not_set_value || vmax != not_set_value) {
    std::vector<T> data;
    filterdata.get(varname, data);
    processWhereMinMax(data, vmin, vmax, where);
  }
}

// -----------------------------------------------------------------------------
template <>
void applyMinMax<util::DateTime>(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                                 ObsFilterData const & filterdata, Variable const & varname) {
  const std::string not_set_value("0000-00-00T00:00:00Z");
  const std::string vmin = mask.getString("minvalue", not_set_value);
  const std::string vmax = mask.getString("maxvalue", not_set_value);

  // Apply mask min/max
  if (vmin != not_set_value || vmax != not_set_value) {
    std::vector<util::DateTime> data;
    filterdata.get(varname, data);
    processWhereMinMax(data, vmin, vmax, where);
  }
}

// -----------------------------------------------------------------------------
/// \brief Process an `any_bit_set_of` keyword in a `where` clause.
///
/// This function sets to `false` all elements of `where` corresponding to elements of `data` in
/// which all bits with indices `bitIndices` are zero. Bits are numbered from 0 starting from the
/// least significant bit.
///
/// The vectors `data` and `where` must be of the same length.
///
/// Example: Suppose `data` is set to [1, 3, 4, 8] and bitIndices to [0, 2]. Then this function will
/// set only the last element of `where` to false, since 8 is the only integer from `data` in whose
/// binary representation both bits 0 and 2 are zero.
void processWhereAnyBitSetOf(const std::vector<int> & data,
                             const std::set<int> & bitIndices,
                             std::vector<bool> & where) {
  std::bitset<32> mask_bs;
  for (const int &bitIndex : bitIndices) {
    mask_bs[bitIndex] = 1;
  }
  const int mask = mask_bs.to_ulong();

  for (size_t jj = 0; jj < data.size(); ++jj) {
    if ((data[jj] & mask) == 0) {
      // None of the specified bits is set
      where[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------
/// \brief Process an `any_bit_unset_of` keyword in a `where` clause.
///
/// This function sets to `false` all elements of `where` corresponding to elements of `data` in
/// which all bits with indices `bitIndices` are non-zero. Bits are numbered from 0 starting from
/// the least significant bit.
///
/// The vectors `data` and `where` must be of the same length.
///
/// Example: Suppose `data` is set to [1, 3, 4, 5] and bitIndices to [0, 2]. Then this function will
/// set only the last element of `where` to false, since 5 is the only integer from `data` in whose
/// binary representation both bits 0 and 2 are non-zero.
void processWhereAnyBitUnsetOf(const std::vector<int> & data,
                               const std::set<int> & bitIndices,
                               std::vector<bool> & where) {
  std::bitset<32> mask_bs;
  for (const int &bitIndex : bitIndices) {
    mask_bs[bitIndex] = 1;
  }
  const int mask = mask_bs.to_ulong();

  for (size_t jj = 0; jj < data.size(); ++jj) {
    if ((data[jj] & mask) == mask) {
      // None of the specified bits is unset
      where[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------
/// \brief Process a `matches_regex` keyword in a `where` clause.
///
/// This function sets to `false` all elements of `where` corresponding to elements of `data` that
/// do not match the regular expression `pattern`. The vectors `data` and `where` must be of the
/// same length.
void processWhereMatchesRegex(const std::vector<std::string> & data,
                              const std::string & pattern,
                              std::vector<bool> & where) {
  std::regex regex(pattern);
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !std::regex_match(data[jj], regex))
      where[jj] = false;
  }
}

/// \brief Process a `matches_regex` keyword in a `where` clause.
///
/// This function sets to `false` all elements of `where` corresponding to elements of `data` whose
/// string representations do not match the regular expression `pattern`. The vectors `data` and
/// `where` must be of the same length.
void processWhereMatchesRegex(const std::vector<int> & data,
                              const std::string & pattern,
                              std::vector<bool> & where) {
  std::regex regex(pattern);
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !std::regex_match(std::to_string(data[jj]), regex))
      where[jj] = false;
  }
}

// -----------------------------------------------------------------------------
/// Returns true if `string` matches any of the patterns from the list `patterns`.
///
/// The patterns may contain wildcards `*` (matching any sequence of characters) and `?` (matching
/// a single character).
bool stringMatchesAnyWildcardPattern(const std::string &string,
                                     const std::vector<std::string> & patterns) {
  return std::any_of(patterns.begin(),
                     patterns.end(),
                     [&string] (const std::string &pattern)
                     { return util::matchesWildcardPattern(string, pattern); });
}

/// \brief Function used to process a `matches_wildcard` or `matches_any_wildcard` keyword in a
/// `where` clause.
///
/// This function sets to `false` all elements of `where` corresponding to elements of `data` that
/// do not match any of the patterns from the list `patterns`. The patterns may contain wildcards
/// `*` (matching any sequence of characters) and `?` (matching a single character). The vectors
/// `data` and `where` must be of the same length.
void processWhereMatchesAnyWildcardPattern(const std::vector<std::string> & data,
                                           const std::vector<std::string> & patterns,
                                           std::vector<bool> & where) {
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !stringMatchesAnyWildcardPattern(data[jj], patterns))
      where[jj] = false;
  }
}

/// \overload Same as the function above, but taking a vector of integers rather than strings.
/// The integers are converted to strings before pattern matching.
void processWhereMatchesAnyWildcardPattern(const std::vector<int> & data,
                                           const std::vector<std::string> & patterns,
                                           std::vector<bool> & where) {
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !stringMatchesAnyWildcardPattern(std::to_string(data[jj]), patterns))
      where[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void isInString(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<std::string> data;
  std::vector<std::string> whitelistvec = mask.getStringVector("is_in");
  std::set<std::string> whitelist(whitelistvec.begin(), whitelistvec.end());
  filterdata.get(varname, data);
  processWhereIsIn(data, whitelist, where);
}

// -----------------------------------------------------------------------------
void isInInteger(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                 ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<int> data;
  std::set<int> whitelist = oops::parseIntSet(mask.getString("is_in"));
  filterdata.get(varname, data);
  processWhereIsIn(data, whitelist, where);
}

// -----------------------------------------------------------------------------
void isNotInString(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                   ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<std::string> data;
  std::vector<std::string> blacklistvec = mask.getStringVector("is_not_in");
  std::set<std::string> blacklist(blacklistvec.begin(), blacklistvec.end());
  filterdata.get(varname, data);
  processWhereIsNotIn(data, blacklist, where);
}

// -----------------------------------------------------------------------------
void isNotInInteger(std::vector<bool> & where, eckit::LocalConfiguration const & mask,
                    ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<int> data;
  filterdata.get(varname, data);
  std::set<int> blacklist = oops::parseIntSet(mask.getString("is_not_in"));
  processWhereIsNotIn(data, blacklist, where);
}

// -----------------------------------------------------------------------------
std::vector<bool> processWhere(const eckit::Configuration & config,
                               const ObsFilterData & filterdata) {
  const size_t nlocs = filterdata.nlocs();

// Everywhere by default if no mask
  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks = config.getSubConfigurations();

  for (size_t jm = 0; jm < masks.size(); ++jm) {
    eckit::LocalConfiguration varconf(masks[jm], "variable");
    Variable var(varconf);
    for (size_t jvar = 0; jvar < var.size(); ++jvar) {
      if (var.group() != "VarMetaData") {
        const Variable varname = var[jvar];
        ioda::ObsDtype dtype = filterdata.dtype(varname);

        if (dtype == ioda::ObsDtype::DateTime) {
          applyMinMax<util::DateTime>(where, masks[jm], filterdata, varname);
        } else if (dtype == ioda::ObsDtype::Integer) {
          applyMinMax<int>(where, masks[jm], filterdata, varname);
        } else {
          applyMinMax<float>(where, masks[jm], filterdata, varname);
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
            isInString(where, masks[jm], filterdata, varname);
          } else if (dtype == ioda::ObsDtype::Integer) {
            isInInteger(where, masks[jm], filterdata, varname);
          } else {
            throw eckit::UserError(
              "Only integer and string variables may be used for processWhere 'is_in'",
              Here());
          }
        }

//      Apply mask is_not_in
        if (masks[jm].has("is_not_in")) {
          if (dtype == ioda::ObsDtype::String) {
            isNotInString(where, masks[jm], filterdata, varname);
          } else if (dtype == ioda::ObsDtype::Integer) {
            isNotInInteger(where, masks[jm], filterdata, varname);
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
            std::set<int> bitIndices = oops::parseIntSet(masks[jm].getString("any_bit_set_of"));
            filterdata.get(varname, data);
            processWhereAnyBitSetOf(data, bitIndices, where);
          } else {
            throw eckit::UserError(
              "Only integer variables may be used for processWhere 'any_bit_set_of'",
              Here());
          }
        }

//      Apply mask any_bit_unset_of
        if (masks[jm].has("any_bit_unset_of")) {
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            std::set<int> bitIndices = oops::parseIntSet(masks[jm].getString("any_bit_unset_of"));
            filterdata.get(varname, data);
            processWhereAnyBitUnsetOf(data, bitIndices, where);
          } else {
            throw eckit::UserError(
              "Only integer variables may be used for processWhere 'any_bit_unset_of'",
              Here());
          }
        }

//      Apply mask matches_regex
        if (masks[jm].has("matches_regex")) {
          const std::string pattern = masks[jm].getString("matches_regex");
          // Select observations for which the variable 'varname' matches the regular expression
          // 'pattern'.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesRegex(data, pattern, where);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesRegex(data, pattern, where);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere 'matches_regex'",
              Here());
          }
        }

//      Apply mask matches_wildcard
        if (masks[jm].has("matches_wildcard")) {
          const std::string pattern = masks[jm].getString("matches_wildcard");
          // Select observations for which the variable 'varname' matches the pattern
          // 'pattern', which may contain the * and ? wildcards.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, {pattern}, where);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, {pattern}, where);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere 'matches_wildcard'",
              Here());
          }
        }

//      Apply mask matches_any_wildcard
        if (masks[jm].has("matches_any_wildcard")) {
          const std::vector<std::string> patterns =
              masks[jm].getStringVector("matches_any_wildcard");
          // Select observations for which the variable 'varname' matches any of the patterns
          // 'patterns'; these may contain the * and ? wildcards.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, patterns, where);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, patterns, where);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere "
              "'matches_any_wildcard'",
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
