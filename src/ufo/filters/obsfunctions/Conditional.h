/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// Parameters for one where configuration from the cases section
/// of the yaml file.
///
/// Example:
///
///   - where:
///     - variable:
///         name: float_variable_2@MetaData
///       minvalue: 0
///     value: 0.5
///
class LocalConditionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalConditionalParameters, Parameters)

 public:
  /// Where clause needed for assignment.  This is passed to create a ProcessWhere
  /// object and requires the expected parameters for ufo::WhereParameters
  oops::RequiredParameter<std::vector<WhereParameters>> where{"where", this};

  /// \brief Value to be assigned when this particular where clause is true.
  oops::RequiredParameter<float> value{"value", this};
};

/// Parameters controlling the Conditional obs function.
class ConditionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ConditionalParameters, Parameters)

 public:
  /// List of cases for assignment.  See LocalConditionalParameters
  /// for what each where clause requires.
  oops::RequiredParameter<std::vector<LocalConditionalParameters>> cases{"cases", this};

  /// Default value for the array to be assigned.  Missing value is used
  /// if this is not present in the yaml.
  oops::OptionalParameter<float> defaultvalue{"defaultvalue", this};

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
/// processWhere class can be used with this obs function.  This can only be used to assign
/// float and int variables because the ObsFunction only computes floating point arrays.
///
/// Example : Create a new integer variable `emissivity@ObsDerived` and assign values based
/// on the surface type.
///
///   - filter: Variable Assignment
///     assignments:
///     - name: emissivity@ObsDerived
///       type: int
///       function:
///         name: Conditional@ObsFunction
///         options:
///           defaultvalue: 0.0 # default value - rttov to calculate.
///           cases:
///           - where:
///             - variable:
///                 name: surface_type@MetaData
///               is_in: 1
///             # if necessary, further conditions could be specified in extra items
///             # in the 'where' list
///             value: 0.3
///           - where:
///             - variable:
///                 name: surface_type@MetaData
///               is_in: 2
///             value: 0.5
///
class Conditional : public ObsFunctionBase {
 public:
  explicit Conditional(const eckit::LocalConfiguration & = eckit::LocalConfiguration());
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ConditionalParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CONDITIONAL_H_
