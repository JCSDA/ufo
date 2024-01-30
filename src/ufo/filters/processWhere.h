/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROCESSWHERE_H_
#define UFO_FILTERS_PROCESSWHERE_H_

#include <set>
#include <string>
#include <vector>

#include "oops/util/AnyOf.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraitsAnyOf.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/PartialDateTime.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
  enum class WhereOperator {AND, OR};

  struct WhereOperatorParameterTraitsHelper {
    typedef WhereOperator EnumType;
    static constexpr char enumTypeName[] = "WhereOperator";
    static constexpr util::NamedEnumerator<WhereOperator> namedValues[] =
      { { WhereOperator::AND, "and" },
        { WhereOperator::OR, "or" }
      };
  };

  enum class WhereValue {VALID, NOT_VALID, TRUE, FALSE};

  struct WhereValueParameterTraitsHelper {
    typedef WhereValue EnumType;
    static constexpr char enumTypeName[] = "WhereValue";
    static constexpr util::NamedEnumerator<WhereValue> namedValues[] =
      { { WhereValue::VALID,     "is_valid" },
        { WhereValue::NOT_VALID, "is_not_valid" },
        { WhereValue::TRUE,      "is_true" },
        { WhereValue::FALSE,     "is_false" }
      };
  };

}  // namespace ufo

namespace oops {

  template<>
  struct ParameterTraits<ufo::WhereOperator> :
    public EnumParameterTraits<ufo::WhereOperatorParameterTraitsHelper>
  {};

  template<>
  struct ParameterTraits<ufo::WhereValue> :
    public EnumParameterTraits<ufo::WhereValueParameterTraitsHelper>
  {};
}  // namespace oops

namespace ufo {
  class ObsFilterData;
  class Variables;

/// \brief Contents of a single element of the list taken by the `where` option of each filter,
/// used to select the subset of observation locations that should be processed by the filter.
class WhereParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WhereParameters, Parameters);

 public:
  /// Variable whose value must fulfil the conditions specified in the remaining parameters.
  oops::RequiredParameter<Variable> variable{"variable", this};

  /// This is a tolerance used for absolute float comparison.  Currently used for float
  /// comparison with isClose and isNotClose.
  oops::OptionalParameter<float> absolutetolerance{"absolute_tolerance", this};

  /// This is a tolerance used for relative tolerance comparison.  Currently used for float
  /// comparison with isClose and isNotClose.
  oops::OptionalParameter<float> relativetolerance{"relative_tolerance", this};

  /// Select locations at which the condition variable is greater than or equal to the specified
  /// value. Can be set to an int, float or datetime in the ISO 8601 format (if any datetime
  /// components are zero, they are ignored).
  oops::OptionalParameter<util::AnyOf<int, float, util::PartialDateTime>> minvalue{
    "minvalue", this};

  /// Select locations at which the condition variable is less than or equal to the specified
  /// value. Can be set to an int, float or datetime in the ISO 8601 format (if any datetime
  /// components are zero, they are ignored).
  oops::OptionalParameter<util::AnyOf<int, float, util::PartialDateTime>> maxvalue{
    "maxvalue", this};

  /// Select whether the minvalue or maxvalue is an exclusive limit. Default value is inclusive.
  oops::Parameter<bool> minExclusive{"min_exclusive", false, this};
  oops::Parameter<bool> maxExclusive{"max_exclusive", false, this};

  /// Select locations at which the condition variable takes one of the specified values. For
  /// integer variables, this can be an integer, range of integers (e.g. `3-5`) or a
  /// comma-separated list of integers and/or ranges (e.g. `3-5, 7-8, 10`). For string variables,
  /// this should be a bracketed list of strings (e.g. `[abc, def]`);
  oops::OptionalParameter<util::AnyOf<std::set<int>, std::vector<std::string>>> isIn{
    "is_in", this};

  /// Select locations at which the condition variable is within tolerance to specified values.
  /// For float variables, this should be a bracketed list of floats (e.g. `[0.0, 0.5]`)
  /// which are compared within a tolerance.  A tolerance must be provided by the user and is a
  /// float either called absolute_tolerance or relative_tolerance (see above);
  oops::OptionalParameter<std::vector<float>> isClose{"is_close_to_any_of", this};

  /// Select locations at which the condition variable does not take any of the specified values.
  /// For integer variables, this can be an integer, range of integers (e.g. `3-5`) or a
  /// comma-separated list of integers and/or ranges (e.g. `3-5, 7-8, 10`). For string variables,
  /// this should be a bracketed list of strings (e.g. `[abc, def]`);
  oops::OptionalParameter<util::AnyOf<std::set<int>, std::vector<std::string>>> isNotIn{
    "is_not_in", this};

  /// Select locations at which the condition variable is not within tolerance of the
  /// specified values. For float variables, this should be a bracketed list of floats
  /// (e.g. `[0.0, 0.5]`) which are compared within a tolerance  A tolerance must be
  /// provided by the user and is a float either called absolute_tolerance or
  /// relative_tolerance (see above);
  oops::OptionalParameter<std::vector<float>> isNotClose{"is_not_close_to_any_of", this};

  /// Select locations at which the condition variable is:
  /// - not set to the missing value indicator or
  /// - set to the missing value indicator or
  /// - is true (typically a diagnostic flag) or
  /// - is false (typically a diagnostic flag)
  oops::OptionalParameter<WhereValue> valueIs{"value", this};

  /// Select locations at which any of the specified bits in the condition variable is set.
  oops::OptionalParameter<std::set<int>> anyBitSetOf{"any_bit_set_of", this};

  /// Select locations at which any of the specified bits in the condition variable is unset.
  oops::OptionalParameter<std::set<int>> anyBitUnsetOf{"any_bit_unset_of", this};

  /// Select locations at which the condition variable matches the specified wildcard (which may
  /// contain the wildcard characters `*`, standing for any number of characters, and `?`, standing
  /// for any character).
  oops::OptionalParameter<std::string> matchesWildcard{"matches_wildcard", this};

  /// Select locations at which the condition variable matches any of the specified wildcards.
  oops::OptionalParameter<std::vector<std::string>> matchesAnyWildcard{
    "matches_any_wildcard", this};

  /// Select locations at which the condition variable matches the specified regular expression.
  oops::OptionalParameter<std::string> matchesRegex{"matches_regex", this};
};

ufo::Variables getAllWhereVariables(const std::vector<WhereParameters> &);
std::vector<bool> processWhere(const std::vector<WhereParameters> &, const ObsFilterData &,
                               const WhereOperator & whereOperator);

}  // namespace ufo

#endif  // UFO_FILTERS_PROCESSWHERE_H_
