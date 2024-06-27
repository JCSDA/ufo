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

#if USE_BOOST_REGEX
#include <boost/regex.hpp>
#define REGEX_NAMESPACE boost
#else
#include <regex>
#define REGEX_NAMESPACE std
#endif

#include "eckit/types/FloatCompare.h"
#include "ioda/ObsSpace.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/wildcard.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

namespace ufo {

constexpr char WhereOperatorParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<WhereOperator>
  WhereOperatorParameterTraitsHelper::namedValues[];

constexpr char WhereValueParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<WhereValue>
  WhereValueParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------
ufo::Variables getAllWhereVariables(const std::vector<WhereParameters> & params) {
  ufo::Variables vars;
  for (const WhereParameters & currentParams : params) {
    vars += currentParams.variable;
  }
  return vars;
}


// -----------------------------------------------------------------------------
template<typename T>
void processWhereMinMax(const std::vector<T> & data,
                        const T & vmin, const T & vmax,
                        std::vector<bool> & mask,
                        bool minExclusive,
                        bool maxExclusive) {
  ASSERT(data.size() == mask.size());
  const T missing = util::missingValue<T>();
  if (vmin != missing || vmax != missing) {
    for (size_t jj = 0; jj < data.size(); ++jj) {
      if (data[jj] == missing) continue;
      if (!minExclusive && vmin != missing && data[jj] < vmin) mask[jj] = false;
      else if (minExclusive && vmin != missing && data[jj] <= vmin) mask[jj] = false;
      if (!maxExclusive && vmax != missing && data[jj] > vmax) mask[jj] = false;
      else if (maxExclusive && vmax != missing && data[jj] >= vmax) mask[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------
void processWhereMinMax(const std::vector<util::DateTime> & data,
                        const util::PartialDateTime & vmin, const util::PartialDateTime & vmax,
                        std::vector<bool> & mask,
                        bool minExclusive,
                        bool maxExclusive) {
  ASSERT(data.size() == mask.size());
  const util::PartialDateTime not_set_value {};
  const util::DateTime missing = util::missingValue<util::DateTime>();
  if (vmin != not_set_value || vmax != not_set_value) {
    for (size_t jj = 0; jj < data.size(); ++jj) {
      if (data[jj] == missing) continue;
      if (!minExclusive && vmin != missing && data[jj] < vmin) mask[jj] = false;
      else if (minExclusive && vmin != missing && data[jj] <= vmin) mask[jj] = false;
      if (!maxExclusive && vmax != missing && data[jj] > vmax) mask[jj] = false;
      else if (maxExclusive && vmax != missing && data[jj] >= vmax) mask[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------
template<typename T>
void processWhereIsDefined(const ObsFilterData & filterdata,
                           const Variable & varname,
                           std::vector<bool> & mask) {
  const T missing = util::missingValue<T>();
  std::vector<T> data;
  filterdata.get(varname, data);
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (data[jj] == missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
template<typename T>
void processWhereIsNotDefined(const ObsFilterData & filterdata,
                              const Variable & varname,
                              std::vector<bool> & mask) {
  const T missing = util::missingValue<T>();
  std::vector<T> data;
  filterdata.get(varname, data);
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (data[jj] != missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
template <class T>
void processWhereIsIn(const std::vector<T> & data,
                      const std::set<T> & whitelist,
                      std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (!oops::contains(whitelist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void processWhereIsClose(const std::vector<float> & data,
                         const float tolerance, const bool relative,
                         const std::vector<float> & whitelist,
                         std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    bool inlist = false;
    for (auto testvalue : whitelist) {
      if (relative) {
        float relativetolerance = testvalue * tolerance;
        if (eckit::types::is_approximately_equal(data[jj], testvalue, relativetolerance)) {
          inlist = true;
          break;
        }
      } else {
        if (eckit::types::is_approximately_equal(data[jj], testvalue, tolerance)) {
          inlist = true;
          break;
        }
      }
    }  // testvalue
    if (!inlist) mask[jj] = false;
  }  // jj
}

// -----------------------------------------------------------------------------
template <class T>
void processWhereIsNotIn(const std::vector<T> & data,
                         const std::set<T> & blacklist,
                         std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  const T missing = util::missingValue<T>();
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (data[jj] == missing || oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void processWhereIsNotIn(const std::vector<std::string> & data,
                         const std::set<std::string> & blacklist,
                         std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void processWhereIsNotClose(const std::vector<float> & data,
                            const float tolerance, const bool relative,
                            const std::vector<float> & blacklist,
                            std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  const float missing = util::missingValue<float>();
  for (size_t jj = 0; jj < data.size(); ++jj) {
    for (auto testvalue : blacklist) {
      if (relative) {
        float relativetolerance = testvalue * tolerance;
        if (data[jj] == missing ||
            eckit::types::is_approximately_equal(data[jj], testvalue, relativetolerance)) {
          mask[jj] = false;
          break;
        }
      } else {
        if (data[jj] == missing ||
            eckit::types::is_approximately_equal(data[jj], testvalue, tolerance)) {
          mask[jj] = false;
          break;
        }
      }
    }  // testvalue
  }  // jj
}

// -----------------------------------------------------------------------------
/// Clear the `mask` wherever `data` is `false`.
void processWhereIsTrue(const std::vector<DiagnosticFlag> & data,
                        std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (!data[jj]) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
/// Clear the `mask` wherever `data` is `true`.
void processWhereIsFalse(const std::vector<DiagnosticFlag> & data,
                         std::vector<bool> & mask) {
  ASSERT(data.size() == mask.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (data[jj]) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
void applyMinMax(std::vector<bool> & where, WhereParameters const & parameters,
                 ObsFilterData const & filterdata, Variable const & varname) {
  const T not_set_value = util::missingValue<T>();

  // Set vmin to the value of the 'minvalue' option if it exists; if not, leave vmin unchanged.
  T vmin = not_set_value;
  if (parameters.minvalue.value() != boost::none)
    vmin = parameters.minvalue.value()->as<T>();
  // Set vmax to the value of the 'maxvalue' option if it exists; if not, leave vmax unchanged.
  T vmax = not_set_value;
  if (parameters.maxvalue.value() != boost::none)
    vmax = parameters.maxvalue.value()->as<T>();

  // Apply mask min/max
  if (vmin != not_set_value || vmax != not_set_value) {
    std::vector<T> data;
    filterdata.get(varname, data);
    processWhereMinMax(data, vmin, vmax, where, parameters.minExclusive, parameters.maxExclusive);
  }
}

// -----------------------------------------------------------------------------
template <>
void applyMinMax<util::DateTime>(std::vector<bool> & where, WhereParameters const & parameters,
                                 ObsFilterData const & filterdata, Variable const & varname) {
  util::PartialDateTime vmin {}, vmax {}, not_set_value {};
  if (parameters.minvalue.value() != boost::none)
    vmin = parameters.minvalue.value()->as<util::PartialDateTime>();
  if (parameters.maxvalue.value() != boost::none)
    vmax = parameters.maxvalue.value()->as<util::PartialDateTime>();

  // Apply mask min/max
  if (vmin != not_set_value || vmax != not_set_value) {
    std::vector<util::DateTime> data;
    filterdata.get(varname, data);
    processWhereMinMax(data, vmin, vmax, where, parameters.minExclusive, parameters.maxExclusive);
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
  ASSERT(data.size() == where.size());
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
  ASSERT(data.size() == where.size());
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
  ASSERT(data.size() == where.size());
  REGEX_NAMESPACE::regex regex(pattern);
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !REGEX_NAMESPACE::regex_match(data[jj], regex))
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
  ASSERT(data.size() == where.size());
  REGEX_NAMESPACE::regex regex(pattern);
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !REGEX_NAMESPACE::regex_match(std::to_string(data[jj]), regex))
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
  ASSERT(data.size() == where.size());
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
  ASSERT(data.size() == where.size());
  for (size_t jj = 0; jj < data.size(); ++jj) {
    if (where[jj] && !stringMatchesAnyWildcardPattern(std::to_string(data[jj]), patterns))
      where[jj] = false;
  }
}

// -----------------------------------------------------------------------------
void isInString(std::vector<bool> & where, std::vector<std::string> const & allowedValues,
                ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<std::string> data;
  std::set<std::string> whitelist(allowedValues.begin(), allowedValues.end());
  filterdata.get(varname, data);
  processWhereIsIn(data, whitelist, where);
}

// -----------------------------------------------------------------------------
void isInInteger(std::vector<bool> & where, std::set<int> const & allowedValues,
                 ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<int> data;
  filterdata.get(varname, data);
  processWhereIsIn(data, allowedValues, where);
}

// -----------------------------------------------------------------------------
void isNotInString(std::vector<bool> & where, std::vector<std::string> const & forbiddenValues,
                   ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<std::string> data;
  std::set<std::string> blacklist(forbiddenValues.begin(), forbiddenValues.end());
  filterdata.get(varname, data);
  processWhereIsNotIn(data, blacklist, where);
}

// -----------------------------------------------------------------------------
void isNotInInteger(std::vector<bool> & where, std::set<int> const & forbiddenValues,
                    ObsFilterData const & filterdata, Variable const & varname) {
  std::vector<int> data;
  filterdata.get(varname, data);
  processWhereIsNotIn(data, forbiddenValues, where);
}

void setWhereVector(std::vector<bool> & where,
                    const bool value) {
  std::fill(where.begin(), where.end(), value);
}

void applyWhereOperator(const WhereOperator & whereOperator,
                        std::vector<bool> & whereTest,
                        std::vector<bool> & where) {
  ASSERT(whereTest.size() == where.size());
  switch (whereOperator) {
  case WhereOperator::AND:
    for (size_t jloc = 0; jloc < whereTest.size(); ++jloc)
      where[jloc] = where[jloc] && whereTest[jloc];
    break;
  case WhereOperator::OR:
    for (size_t jloc = 0; jloc < whereTest.size(); ++jloc)
      where[jloc] = where[jloc] || whereTest[jloc];
    break;
  }
}

// -----------------------------------------------------------------------------
std::vector<bool> processWhere(const std::vector<WhereParameters> & params,
                               const ObsFilterData & filterdata,
                               const WhereOperator & whereOperator) {
  const size_t nlocs = filterdata.nlocs();

  // Vector to which all selection operations are applied.
  std::vector<bool> where(nlocs);
  switch (whereOperator) {
  case WhereOperator::AND:
    // Set `where` to `true` everywhere because a logical `and` will be used.
    setWhereVector(where, true);
    break;
  case WhereOperator::OR:
    // Set `where` to `false` everywhere because a logical `or` will be used.
    setWhereVector(where, false);
    break;
  }

  // Set `where` to `true` everywhere if there are no comparisons to be made.
  if (params.empty())
    setWhereVector(where, true);

  for (const WhereParameters &currentParams : params) {
    const Variable &var = currentParams.variable;
    for (size_t jvar = 0; jvar < var.size(); ++jvar) {
        const Variable varname = var[jvar];
        ioda::ObsDtype dtype = filterdata.dtype(varname);

        // Vector to which each operation is applied individually.
        // Set to `true` at the start of each test.
        std::vector<bool> whereTest(nlocs);

//      Apply mask min/max
        if (currentParams.minvalue.value() ||
            currentParams.maxvalue.value()) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::DateTime) {
            applyMinMax<util::DateTime>(whereTest, currentParams, filterdata, varname);
          } else if (dtype == ioda::ObsDtype::Integer) {
            applyMinMax<int>(whereTest, currentParams, filterdata, varname);
          } else {
            applyMinMax<float>(whereTest, currentParams, filterdata, varname);
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_defined
        if (currentParams.valueIs.value() &&
            (*currentParams.valueIs.value() == WhereValue::VALID)) {
          setWhereVector(whereTest, true);
          if (filterdata.has(varname)) {
            if (dtype == ioda::ObsDtype::Integer) {
              processWhereIsDefined<int>(filterdata, varname, whereTest);
            } else if (dtype == ioda::ObsDtype::Float) {
              processWhereIsDefined<float>(filterdata, varname, whereTest);
            } else if (dtype == ioda::ObsDtype::String) {
              processWhereIsDefined<std::string>(filterdata, varname, whereTest);
            } else if (dtype != ioda::ObsDtype::Empty) {
              throw eckit::UserError(
                "Only integer, float and string variables may be used for processWhere "
                "'is_defined'",
                Here());
            }
          } else {
            std::fill(whereTest.begin(), whereTest.end(), false);
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_not_defined
        if (currentParams.valueIs.value() &&
            (*currentParams.valueIs.value() == WhereValue::NOT_VALID)) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::Integer) {
            processWhereIsNotDefined<int>(filterdata, varname, whereTest);
          } else if (dtype == ioda::ObsDtype::Float) {
            processWhereIsNotDefined<float>(filterdata, varname, whereTest);
          } else if (dtype == ioda::ObsDtype::String) {
            processWhereIsNotDefined<std::string>(filterdata, varname, whereTest);
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only integer, float and string variables may be used for processWhere "
              "'is_not_defined'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_in
        if (currentParams.isIn.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::String) {
            isInString(whereTest, currentParams.isIn.value()->as<std::vector<std::string>>(),
                       filterdata, varname);
          } else if (dtype == ioda::ObsDtype::Integer) {
            isInInteger(whereTest, currentParams.isIn.value()->as<std::set<int>>(),
                        filterdata, varname);
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only integer and string variables may be used for processWhere 'is_in'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_close
        if (currentParams.isClose.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::Float) {
            std::vector<float> data;
            filterdata.get(varname, data);
            if (currentParams.relativetolerance.value() == boost::none &&
                currentParams.absolutetolerance.value() != boost::none) {
              processWhereIsClose(data, currentParams.absolutetolerance.value().get(),
                                  false, currentParams.isClose.value().get(), whereTest);
            } else if (currentParams.relativetolerance.value() != boost::none &&
                       currentParams.absolutetolerance.value() == boost::none) {
              processWhereIsClose(data, currentParams.relativetolerance.value().get(),
                                  true, currentParams.isClose.value().get(), whereTest);
            } else {
              throw eckit::UserError(
                "For 'is_close' one (and only one) tolerance is needed.",
                Here());
            }
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only float variables may be used for processWhere 'is_close'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_not_in
        if (currentParams.isNotIn.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::String) {
            isNotInString(whereTest, currentParams.isNotIn.value()->as<std::vector<std::string>>(),
                          filterdata, varname);
          } else if (dtype == ioda::ObsDtype::Integer) {
            isNotInInteger(whereTest, currentParams.isNotIn.value()->as<std::set<int>>(),
                           filterdata, varname);
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only integer and string variables may be used for processWhere 'is_not_in'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_not_close
        if (currentParams.isNotClose.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::Float) {
            std::vector<float> data;
            filterdata.get(varname, data);
            if (currentParams.relativetolerance.value() == boost::none &&
                currentParams.absolutetolerance.value() != boost::none) {
              processWhereIsNotClose(data, currentParams.absolutetolerance.value().get(),
                                     false, currentParams.isNotClose.value().get(), whereTest);
            } else if (currentParams.relativetolerance.value() != boost::none &&
                       currentParams.absolutetolerance.value() == boost::none) {
              processWhereIsNotClose(data, currentParams.relativetolerance.value().get(),
                                     true, currentParams.isNotClose.value().get(), whereTest);
            } else {
              throw eckit::UserError(
                "For 'is_close' one (and only one) tolerance is needed.",
                Here());
            }
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only float variables may be used for processWhere 'is_not_close'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_set
        if (currentParams.valueIs.value() &&
            (*currentParams.valueIs.value() == WhereValue::TRUE)) {
          setWhereVector(whereTest, true);
          if (filterdata.has(varname)) {
            std::vector<DiagnosticFlag> data;
            filterdata.get(varname, data);
            processWhereIsTrue(data, whereTest);
          } else {
            std::fill(whereTest.begin(), whereTest.end(), false);
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask is_not_set
        if (currentParams.valueIs.value() &&
            (*currentParams.valueIs.value() == WhereValue::FALSE)) {
          setWhereVector(whereTest, true);
          std::vector<DiagnosticFlag> data;
          filterdata.get(varname, data);
          processWhereIsFalse(data, whereTest);
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask any_bit_set_of
        if (currentParams.anyBitSetOf.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            const std::set<int> &bitIndices = *currentParams.anyBitSetOf.value();
            filterdata.get(varname, data);
            processWhereAnyBitSetOf(data, bitIndices, whereTest);
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only integer variables may be used for processWhere 'any_bit_set_of'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask any_bit_unset_of
        if (currentParams.anyBitUnsetOf.value() != boost::none) {
          setWhereVector(whereTest, true);
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            const std::set<int> &bitIndices = *currentParams.anyBitUnsetOf.value();
            filterdata.get(varname, data);
            processWhereAnyBitUnsetOf(data, bitIndices, whereTest);
          } else if (dtype != ioda::ObsDtype::Empty) {
            throw eckit::UserError(
              "Only integer variables may be used for processWhere 'any_bit_unset_of'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask matches_regex
        if (currentParams.matchesRegex.value() != boost::none) {
          setWhereVector(whereTest, true);
          const std::string pattern = *currentParams.matchesRegex.value();
          // Select observations for which the variable 'varname' matches the regular expression
          // 'pattern'.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesRegex(data, pattern, whereTest);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesRegex(data, pattern, whereTest);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere 'matches_regex'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask matches_wildcard
        if (currentParams.matchesWildcard.value() != boost::none) {
          setWhereVector(whereTest, true);
          const std::string &pattern = *currentParams.matchesWildcard.value();
          // Select observations for which the variable 'varname' matches the pattern
          // 'pattern', which may contain the * and ? wildcards.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, {pattern}, whereTest);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, {pattern}, whereTest);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere 'matches_wildcard'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
        }

//      Apply mask matches_any_wildcard
        if (currentParams.matchesAnyWildcard.value() != boost::none) {
          setWhereVector(whereTest, true);
          const std::vector<std::string> &patterns = *currentParams.matchesAnyWildcard.value();
          // Select observations for which the variable 'varname' matches any of the patterns
          // 'patterns'; these may contain the * and ? wildcards.
          if (dtype == ioda::ObsDtype::Integer) {
            std::vector<int> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, patterns, whereTest);
          } else if (dtype == ioda::ObsDtype::String) {
            std::vector<std::string> data;
            filterdata.get(varname, data);
            processWhereMatchesAnyWildcardPattern(data, patterns, whereTest);
          } else {
            throw eckit::UserError(
              "Only string and integer variables may be used for processWhere "
              "'matches_any_wildcard'",
              Here());
          }
          applyWhereOperator(whereOperator, whereTest, where);
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
