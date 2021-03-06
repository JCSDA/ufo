/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/VariableAssignment.h"

#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/processWhere.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

namespace {

/// Convert \p valueAsString to `VariableType` and for each vector in \p values assign that value
/// at locations selected by the `where` statement.
template <typename VariableType>
void assignValue(const std::string &valueAsString,
                 const std::vector<bool> &apply,
                 ioda::ObsDataVector<VariableType> &values) {
  VariableType newValue;
  if (!boost::conversion::try_lexical_convert(valueAsString, newValue))
    throw eckit::BadCast("Value '" + valueAsString +
                         "' could not be converted to the required type", Here());

  for (size_t ival = 0; ival < values.nvars(); ++ival) {
    std::vector<VariableType> &currentValues = values[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc)
      if (apply[iloc])
        currentValues[iloc] = newValue;
  }
}

/// Evaluate the ObsFunction \p function and assign the vectors it produced to successive vectors
/// in \p values (only at locations selected by the `where` statement).
template <typename VariableType>
void assignFunction(const ufo::Variable &function,
                    const ufo::Variable &variable,
                    const std::vector<bool> &apply,
                    const ObsFilterData &data,
                    ioda::ObsDataVector<VariableType> &values) {
  ioda::ObsDataVector<float> newValues(data.obsspace(), variable.toOopsVariables());
  data.get(function, newValues);

  const VariableType missing = util::missingValue(VariableType());
  const float missingfloat = util::missingValue(float());
  for (size_t ival = 0; ival < values.nvars(); ++ival) {
    std::vector<VariableType> &currentValues = values[ival];
    const ioda::ObsDataRow<float> &currentNewValues = newValues[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc) {
      if (apply[iloc] && currentNewValues[iloc] != missingfloat)
        currentValues[iloc] = static_cast<VariableType>(currentNewValues[iloc]);
      if (apply[iloc] && currentNewValues[iloc] == missingfloat)
        currentValues[iloc] = missing;
    }
  }
}

/// Assign values to a numeric variable (of type float or int).
template <typename VariableType>
void assignNumericValues(const AssignmentParameters &params,
                         const ufo::Variable &variable,
                         const std::vector<bool> &apply,
                         const ObsFilterData &data,
                         ioda::ObsDataVector<VariableType> &values) {
  if (params.value_.value() != boost::none) {
    assignValue(*params.value_.value(), apply, values);
  } else {
    ASSERT(params.function.value() != boost::none);
    assignFunction(*params.function.value(), variable, apply, data, values);
  }
}

/// Assign values to a non-numeric variable (of type string or DateTime).
template <typename VariableType>
void assignNonnumericValues(const AssignmentParameters &params,
                            const std::vector<bool> &apply,
                            ioda::ObsDataVector<VariableType> &values) {
  if (params.value_.value() != boost::none) {
    assignValue(*params.value_.value(), apply, values);
  } else {
    ASSERT(params.function.value() != boost::none);
    throw eckit::BadValue("ObsFunction values cannot be assigned to non-numeric variables", Here());
  }
}

/// Retrieve and return the current values of the variable \p variable from \p obsdb (as
/// vectors). Variables that aren't currently stored in \p obsdb are treated as if they consisted
/// entirely of missing values.
template <typename VariableType>
ioda::ObsDataVector<VariableType> getCurrentValues(const ufo::Variable &variable,
                                                   ioda::ObsSpace &obsdb) {
  ioda::ObsDataVector<VariableType> values(obsdb, variable.toOopsVariables());
  for (size_t ich = 0; ich < variable.size(); ++ich) {
    const std::string variableWithChannel = variable.variable(ich);
    if (obsdb.has(variable.group(), variableWithChannel)) {
      // Variable exists -- retrieve its values from the ObsSpace
      obsdb.get_db(variable.group(), variableWithChannel, values[ich]);
    } else {
      // Variable doesn't exist yet -- fill the vector with missing values
      values[ich].assign(obsdb.nlocs(), util::missingValue(VariableType()));
    }
  }
  return values;
}

/// Save the values \p values of variable \p variable to \p obsdb.
template <typename VariableType>
void saveValues(const ufo::Variable &variable,
                const ioda::ObsDataVector<VariableType> &values,
                ioda::ObsSpace &obsdb) {
  for (size_t ich = 0; ich < variable.size(); ++ich)
    obsdb.put_db(variable.group(), variable.variable(ich), values[ich]);
}

/// Retrieve the current values of a numeric variable \p variable from \p obsdb (or if it doesn't
/// already exist, fill it with missing values), assign new values to elements selected by the
/// `where` clause and save the results to \p obsdb.
template <typename VariableType>
void assignToNumericVariable(const ufo::Variable &variable,
                             const AssignmentParameters &params,
                             const std::vector<bool> &apply,
                             const ObsFilterData &data,
                             ioda::ObsSpace &obsdb) {
  ioda::ObsDataVector<VariableType> values = getCurrentValues<VariableType>(variable, obsdb);
  assignNumericValues(params, variable, apply, data, values);
  saveValues(variable, values, obsdb);
}

/// Retrieve the current values of a non-numeric variable \p variable from \p obsdb (or if it
/// doesn't already exist, fill it with missing values), assign new values to elements selected by
/// the `where` clause and save the results to \p obsdb.
template <typename VariableType>
void assignToNonnumericVariable(const ufo::Variable &variable,
                                const AssignmentParameters &params,
                                const std::vector<bool> &apply,
                                ioda::ObsSpace &obsdb) {
  ioda::ObsDataVector<VariableType> values = getCurrentValues<VariableType>(variable, obsdb);
  assignNonnumericValues(params, apply, values);
  saveValues(variable, values, obsdb);
}

/// Delegate work to an appropriate function depending on whether \p dtype is a numeric
/// or non-numeric type.
void assignToVariable(const ufo::Variable &variable,
                      ioda::ObsDtype dtype,
                      const AssignmentParameters &params,
                      const std::vector<bool> &apply,
                      const ObsFilterData &data,
                      ioda::ObsSpace &obsdb) {
  switch (dtype) {
  case ioda::ObsDtype::Float:
    assignToNumericVariable<float>(variable, params, apply, data, obsdb);
    break;
  case ioda::ObsDtype::Integer:
    assignToNumericVariable<int>(variable, params, apply, data, obsdb);
    break;
  case ioda::ObsDtype::String:
    assignToNonnumericVariable<std::string>(variable, params, apply, obsdb);
    break;
  case ioda::ObsDtype::DateTime:
    assignToNonnumericVariable<util::DateTime>(variable, params, apply, obsdb);
    break;
  default:
    ASSERT_MSG(false, "Unrecognized data type");
  }
}

/// Return the variable to which new values will be assigned.
ufo::Variable getVariable(const AssignmentParameters &params) {
  const std::set<int> setChannels = oops::parseIntSet(params.channels);
  std::vector<int> vecChannels(setChannels.begin(), setChannels.end());
  const ufo::Variable variable(params.name, vecChannels);
  if (variable.group() == "ObsValue")
    throw eckit::BadValue("Assignment to variables from the ObsValue group is not allowed",
                          Here());
  return variable;
}

/// Return the data type of the variable to which new values will be assigned.
ioda::ObsDtype getDataType(boost::optional<ioda::ObsDtype> dtypeParam,
                           const ufo::Variable &variable,
                           const ioda::ObsSpace &obsdb) {
  if (dtypeParam != boost::none) {
    // If the dtype option has been set, return its value.
    return *dtypeParam;
  } else {
    // Otherwise, check if the variable to which new values should be assigned already
    // exists and if so, return its data type.
    for (size_t ich = 0; ich < variable.size(); ++ich) {
      const std::string variableWithChannel = variable.variable(ich);
      if (obsdb.has(variable.group(), variableWithChannel))
        return obsdb.dtype(variable.group(), variableWithChannel);
    }
    // The variable doesn't exist yet.
    throw eckit::BadParameter("You need to specify the type of the variable to be created "
                              "by setting the 'type' option of the filter to 'float', 'int', "
                              "'string' or 'datetime'.");
  }
}

}  // namespace


void AssignmentParameters::deserialize(util::CompositePath &path,
                                       const eckit::Configuration &config) {
  oops::Parameters::deserialize(path, config);

  // These checks should really be done at the validation stage (using JSON Schema),
  // but this isn't supported yet, so this is better than nothing.
  if ((value_.value() == boost::none && function.value() == boost::none) ||
      (value_.value() != boost::none && function.value() != boost::none))
    throw eckit::UserError(path.path() +
                           ": Exactly one of the 'value' and 'function' options must be present");
}


VariableAssignment::VariableAssignment(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : ObsProcessorBase(obsdb, parameters.deferToPost, std::move(flags), std::move(obserr)),
    parameters_(parameters)
{
  oops::Log::debug() << "VariableAssignment: config = " << parameters_ << std::endl;
  allvars_ += getAllWhereVariables(parameters.where);

  for (const AssignmentParameters &assignment : parameters.assignments.value()) {
    if (assignment.function.value() != boost::none) {
      allvars_ += *assignment.function.value();
    }
  }
}

void VariableAssignment::doFilter() const {
  oops::Log::trace() << "VariableAssignment doFilter begin" << std::endl;

  // Select locations at which the filter will be applied
  const std::vector<bool> apply = processWhere(parameters_.where, data_);

  // Assign values to successive sets of variables
  for (const AssignmentParameters &assignment : parameters_.assignments.value()) {
    const ufo::Variable variable = getVariable(assignment);
    const ioda::ObsDtype dtype = getDataType(assignment.type, variable, obsdb_);
    assignToVariable(variable, dtype, assignment, apply, data_, obsdb_);
  }

  oops::Log::trace() << "VariableAssignment doFilter end" << std::endl;
}

void VariableAssignment::print(std::ostream & os) const {
  os << "VariableAssignment: config = " << parameters_ << std::endl;
}

}  // namespace ufo
