/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

/// Parameters for one where configuration from the cases section
/// of the yaml file.
///
/// Example:
///
///   - where:
///     - variable:
///         name: MetaData/float_variable_2
///       minvalue: 0
///     value: 0.5
///
template <typename FunctionValue>
class LocalConditionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalConditionalParameters, Parameters)

 public:
  /// Where clause needed for assignment.  This is passed to create a ProcessWhere
  /// object and requires the expected parameters for ufo::WhereParameters
  oops::RequiredParameter<std::vector<WhereParameters>> where{"where", this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// \brief Value to be assigned when this particular where clause is true.
  oops::RequiredParameter<FunctionValue> value{"value", this};
};

/// Parameters controlling the Conditional obs function.
template <typename FunctionValue>
class ConditionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ConditionalParameters, Parameters)

 public:
  /// List of cases for assignment.  See LocalConditionalParameters
  /// for what each where clause requires.
  oops::RequiredParameter<std::vector<LocalConditionalParameters<FunctionValue>>> cases{
      "cases", this};

  /// Default value for the array to be assigned.  Missing value is used
  /// if this is not present in the yaml.
  oops::OptionalParameter<FunctionValue> defaultvalue{"defaultvalue", this};

  /// When this flag is true a value is assigned for the first matching where case
  /// for a location.  This matches with the python case logic.
  /// When this is false the last matching case will take precedence which
  /// replicates a series of separate variable assignment filter calls.
  oops::Parameter<bool> firstmatchingcase{"firstmatchingcase", true, this};
};

/// \brief Creates an array with values for specified variables selected by a series of where
/// statements.
///
/// The obs function has been designed primarily to work with the Variable assignment filter
/// to simplify the assignment of more complicated variables.  Any functionality in the
/// processWhere class can be used with this obs function.
///
/// This template is used to define four ObsFunctions, each producing values of a different type:
/// * `ObsFunction/Conditional` produces floats
/// * `IntObsFunction/Conditional` produces ints
/// * `StringObsFunction/Conditional` produces strings
/// * `DateTimeObsFunction/Conditional` produces datetimes.
///
/// Example 1: Create a new floating-point variable `ObsDerived/emissivity` and assign values based
/// on the surface type.
///
///     - filter: Variable Assignment
///       assignments:
///       - name: ObsDerived/emissivity
///         type: float
///         function:
///           name: ObsFunction/Conditional
///           options:
///             defaultvalue: 0.0 # default value - rttov to calculate.
///             cases:
///             - where:
///               - variable:
///                   name: MetaData/surfaceQualifier
///                 is_in: 1
///               # if necessary, further conditions could be specified in extra items
///               # in the 'where' list
///               value: 0.3
///             - where:
///               - variable:
///                   name: MetaData/surfaceQualifier
///                 is_in: 2
///               value: 0.5
///
/// Example 2: Create a new string variable `MetaData/surfaceDescription` and set it to `land`,
/// `sea` or `unknown` depending on the value of the `MetaData/surfaceQualifier` variable.
///
///     - filter: Variable Assignment
///       assignments:
///       - name: MetaData/surfaceDescription
///         type: string
///         function:
///           name: StringObsFunction/Conditional
///           options:
///             defaultvalue: unknown
///             cases:
///             - where:
///               - variable:
///                   name: MetaData/surfaceQualifier
///                 is_in: 1
///               value: land
///             - where:
///               - variable:
///                   name: MetaData/surfaceQualifier
///                 is_in: 2
///               value: sea

///
template <typename FunctionValue>
class Conditional : public ObsFunctionBase<FunctionValue> {
 public:
  explicit Conditional(const eckit::LocalConfiguration & = eckit::LocalConfiguration());
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ConditionalParameters<FunctionValue> options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_
