/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/VariableAssignment.h"

#include <algorithm>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

namespace {

/// Convert a float \p x to an int by rounding. If \p is equal to \p missingIn, return \p
/// missingOut. If the value to be returned is too large to be represented by an int, throw an
/// exception.
int safeCast(float x, float missingIn, int missingOut) {
  if (x == missingIn) {
    return missingOut;
  }
  // Throws boost::math::rounding_error if x is outside the range representable by ints
  return boost::math::iround(x);
}

/// Cast an int \p x to a float. If \p is equal to \p missingIn, return \p missingOut.
float safeCast(int x, int missingIn, float missingOut) {
  if (x == missingIn) {
    return missingOut;
  }
  return x;
}

/// Convert a util::DateTime \p x to a numeric type DestinationVariableType by computing the number
/// of seconds relative to a chosen epoch. If \p x is equal to \p missingIn, return \p missingOut.
/// If the number of seconds is too large to be represented by DestinationVariableType,
/// throw an exception. A more appropriate epoch should be chosen in such cases.
template <typename VariableType>
VariableType secondsFromEpoch(const util::DateTime & x,
                              const util::DateTime & missingIn,
                              VariableType missingOut,
                              const util::DateTime & epoch,
                              int64_t minVariableType,
                              int64_t maxVariableType) {
  if (x == missingIn) {
    return missingOut;
  }
  const int64_t nSeconds = (x - epoch).toSeconds();
  if (nSeconds < minVariableType || nSeconds > maxVariableType) {
    throw eckit::BadCast("Invalid epoch causing overflow.", Here());
  }
  return static_cast<VariableType>(nSeconds);
}

/// Convert \p valueAsString to `VariableType` and for each vector in \p values assign that value
/// at locations selected by the `where` statement.
template <typename VariableType>
void assignValue(const std::string &valueAsString,
                 const std::vector<bool> &apply,
                 ioda::ObsDataVector<VariableType> &values) {
  VariableType newValue;
  if (valueAsString == "missing") {
    newValue = util::missingValue<VariableType>();
  } else {
    if (!boost::conversion::try_lexical_convert(valueAsString, newValue)) {
      throw eckit::BadCast("Value '" + valueAsString +
                           "' could not be converted to the required type", Here());
    }
  }

  for (size_t ival = 0; ival < values.nvars(); ++ival) {
    std::vector<VariableType> &currentValues = values[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc)
      if (apply[iloc]) {
        currentValues[iloc] = newValue;
      }
  }
}

/// For each location selected by the `where` statement, copy the corresponding element of
/// \p source to \p destination.
///
/// This two-parameter template is selected by the compiler when SourceVariableType is
/// different from DestinationVariableType. It takes care of converting missing value indicators
/// of type SourceVariableType to ones of type DestinationVariableType.
template <typename SourceVariableType, typename DestinationVariableType>
void assignObsDataVector(const std::vector<bool> &apply,
                         const ioda::ObsDataVector<SourceVariableType> &source,
                         ioda::ObsDataVector<DestinationVariableType> &destination) {
  ASSERT(source.nvars() == destination.nvars());

  const SourceVariableType missingSource = util::missingValue<SourceVariableType>();
  const DestinationVariableType missingDestination = util::missingValue<DestinationVariableType>();

  for (size_t ival = 0; ival < source.nvars(); ++ival) {
    const ioda::ObsDataRow<SourceVariableType> &currentSource = source[ival];
    ioda::ObsDataRow<DestinationVariableType> &currentDestination = destination[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc) {
      if (apply[iloc]) {
        currentDestination[iloc] = safeCast(currentSource[iloc], missingSource, missingDestination);
      }
    }
  }
}

/// For each location selected by the `where` statement, copy the corresponding element of
/// \p source to \p destination. Used when converting util::DateTimes to numerical values.
template <typename VariableType>
void assignObsDataVector(const std::vector<bool> &apply,
                         const ioda::ObsDataVector<util::DateTime> &source,
                         ioda::ObsDataVector<VariableType> &destination,
                         const util::DateTime & epoch) {
  ASSERT(source.nvars() == destination.nvars());

  const util::DateTime missingSource = util::missingValue<util::DateTime>();
  const VariableType missingDestination = util::missingValue<VariableType>();
  // Minimum and maximum destination type cast to int64_t.
  // This ensures that out-of-bounds values (caused by an inappropriate epoch)
  // are not silenty stored in the output vector.
  const int64_t minVariableType =
    static_cast<int64_t>(std::numeric_limits<VariableType>::lowest());
  const int64_t maxVariableType =
    static_cast<int64_t>(std::numeric_limits<VariableType>::max());
  for (size_t ival = 0; ival < source.nvars(); ++ival) {
    const ioda::ObsDataRow<util::DateTime> &currentSource = source[ival];
    ioda::ObsDataRow<VariableType> &currentDestination = destination[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc) {
      if (apply[iloc]) {
        currentDestination[iloc] = secondsFromEpoch(currentSource[iloc],
                                                    missingSource,
                                                    missingDestination,
                                                    epoch,
                                                    minVariableType,
                                                    maxVariableType);
      }
    }
  }
}

/// For each location selected by the `where` statement, copy the corresponding element of
/// \p source to \p destination.
///
/// This single-parameter template is selected by the compiler when \p source and \p destination
/// are of the same type.
template <typename VariableType>
void assignObsDataVector(const std::vector<bool> &apply,
                         const ioda::ObsDataVector<VariableType> &source,
                         ioda::ObsDataVector<VariableType> &destination) {
  ASSERT(source.nvars() == destination.nvars());

  for (size_t ival = 0; ival < source.nvars(); ++ival) {
    const ioda::ObsDataRow<VariableType> &currentSource = source[ival];
    ioda::ObsDataRow<VariableType> &currentDestination = destination[ival];
    for (size_t iloc = 0; iloc < apply.size(); ++iloc) {
      if (apply[iloc]) {
        currentDestination[iloc] = currentSource[iloc];
      }
    }
  }
}

/// Retrieve the variable \p variable of type \c SourceVariableType and assign its components
/// (channels) to successive vectors in \p values (only at locations selected by the `where`
/// statement).
template <typename SourceVariableType, typename DestinationVariableType>
void assignVariable(const ufo::Variable &variable,
                    const bool skipDerived,
                    const std::vector<bool> &apply,
                    const ObsFilterData &data,
                    ioda::ObsDataVector<DestinationVariableType> &values) {
  ioda::ObsDataVector<SourceVariableType> newValues(data.obsspace(), variable.toOopsObsVariables());
  data.get(variable, newValues, skipDerived);
  assignObsDataVector(apply, newValues, values);
}

/// Assign util::DateTimes to numerical values relative to an epoch.
template <typename VariableType>
void assignVariable(const ufo::Variable &variable,
                    const util::DateTime &epoch,
                    const bool skipDerived,
                    const std::vector<bool> &apply,
                    const ObsFilterData &data,
                    ioda::ObsDataVector<VariableType> &values) {
  ioda::ObsDataVector<util::DateTime> newValues(data.obsspace(), variable.toOopsObsVariables());
  data.get(variable, newValues, skipDerived);
  assignObsDataVector(apply, newValues, values, epoch);
}

/// Evaluate the ObsFunction \p function and assign the vectors it produced to successive vectors
/// in \p values (only at locations selected by the `where` statement).
template <typename FunctionValueType, typename VariableType>
void assignFunction(const ufo::Variable &function,
                    const ufo::Variable &variable,
                    const std::vector<bool> &apply,
                    const ObsFilterData &data,
                    ioda::ObsDataVector<VariableType> &values) {
  ioda::ObsDataVector<FunctionValueType> newValues(data.obsspace(), variable.toOopsObsVariables());
  data.get(function, newValues);
  assignObsDataVector(apply, newValues, values);
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
  } else if (params.sourceVariable.value() != boost::none) {
    switch (data.dtype(*params.sourceVariable.value())) {
    case ioda::ObsDtype::Float:
      assignVariable<float>(*params.sourceVariable.value(), params.skipDerived,
                            apply, data, values);
      break;
    case ioda::ObsDtype::Integer:
      assignVariable<int>(*params.sourceVariable.value(), params.skipDerived,
                          apply, data, values);
      break;
    case ioda::ObsDtype::DateTime:
      if (params.epoch.value() != boost::none) {
        assignVariable(*params.sourceVariable.value(),
                       *params.epoch.value(),
                       params.skipDerived,
                       apply, data, values);
      } else {
        throw eckit::UserError("Converting a DateTime to a numeric value requires the "
                               "`epoch` parameter to be set", Here());
      }
      break;
    case ioda::ObsDtype::Empty:
      oops::Log::info() << "ufo::VariableAssignment::assignNumericValues "
                        << "not performed on empty MPI ranks " << std::endl;
      break;
    default:
      throw eckit::BadParameter(params.sourceVariable.value()->fullName() +
                                " is not a numeric variable", Here());
    }
  } else {
    ASSERT(params.function.value() != boost::none);
    if (params.function.value()->group() == ObsFunctionTraits<float>::groupName) {
      assignFunction<float>(*params.function.value(), variable, apply, data, values);
    } else if (params.function.value()->group() == ObsFunctionTraits<int>::groupName) {
      assignFunction<int>(*params.function.value(), variable, apply, data, values);
    } else {
      throw eckit::BadParameter(params.function.value()->fullName() +
                                " is not a function producing numeric values", Here());
    }
  }
}

/// Assign values to a non-numeric variable (of type string or DateTime).
template <typename VariableType>
void assignNonnumericValues(const AssignmentParameters &params,
                            const ufo::Variable &variable,
                            const std::vector<bool> &apply,
                            const ObsFilterData &data,
                            ioda::ObsDataVector<VariableType> &values) {
  if (params.value_.value() != boost::none) {
    assignValue(*params.value_.value(), apply, values);
  } else if (params.sourceVariable.value() != boost::none) {
    assignVariable<VariableType>(*params.sourceVariable.value(), params.skipDerived,
                                 apply, data, values);
  } else {
    ASSERT(params.function.value() != boost::none);
    assignFunction<VariableType>(*params.function.value(), variable, apply, data, values);
  }
}

/// Retrieve and return the current values of the variable \p variable from \p obsdb (as
/// vectors). Variables that aren't currently stored in \p obsdb are treated as if they consisted
/// entirely of missing values.
template <typename VariableType>
ioda::ObsDataVector<VariableType> getCurrentValues(const ufo::Variable &variable,
                                                   ioda::ObsSpace &obsdb,
                                                   const bool skipDerived) {
  ioda::ObsDataVector<VariableType> values(obsdb, variable.toOopsObsVariables());
  for (size_t ich = 0; ich < variable.size(); ++ich) {
    const std::string variableWithChannel = variable.variable(ich);
    if (obsdb.has(variable.group(), variableWithChannel)) {
      // Variable exists -- retrieve its values from the ObsSpace
      obsdb.get_db(variable.group(), variableWithChannel, values[ich], {}, skipDerived);
    } else {
      // Variable doesn't exist yet -- fill the vector with missing values
      values[ich].assign(obsdb.nlocs(), util::missingValue<VariableType>());
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

/// Change the QC flag from `miss` to `pass` if the obs value is no longer missing or from `pass` to
/// `miss` if the obs value is now missing.
void updateQCFlags(const ioda::ObsDataVector<float> &obsvalues, ioda::ObsDataVector<int> &qcflags) {
  const float missing = util::missingValue<float>();

  for (size_t ivar = 0; ivar < obsvalues.nvars(); ++ivar) {
    if (qcflags.varnames().has(obsvalues.varnames()[ivar])) {
      const ioda::ObsDataRow<float> &currentValues = obsvalues[ivar];
      ioda::ObsDataRow<int> &currentFlags = qcflags[obsvalues.varnames()[ivar]];
      for (size_t iloc = 0; iloc < obsvalues.nlocs(); ++iloc) {
        if (currentFlags[iloc] == QCflags::missing && currentValues[iloc] != missing) {
          currentFlags[iloc] = QCflags::pass;
        } else if (currentFlags[iloc] == QCflags::pass && currentValues[iloc] == missing) {
          currentFlags[iloc] = QCflags::missing;
        }
      }
    }
  }
}

/// Retrieve the current values of an int-valued variable \p variable from \p obsdb (or if it
/// doesn't already exist, fill it with missing values), assign new values to elements selected by
/// the `where` clause and save the results to \p obsdb.
void assignToIntVariable(const ufo::Variable &variable,
                         const AssignmentParameters &params,
                         const std::vector<bool> &apply,
                         const ObsFilterData &data,
                         ioda::ObsSpace &obsdb) {
  ioda::ObsDataVector<int> values = getCurrentValues<int>(variable, obsdb, params.skipDerived);
  assignNumericValues(params, variable, apply, data, values);
  saveValues(variable, values, obsdb);
}

/// Works like `assignToIntVariable()`, but in addition if \p variable belongs to the `ObsValue` or
/// `DerivedObsValue` group and is one of the simulated variables, it updates `qcflags` by (a)
/// changing `miss` to `pass` if the obs value is no longer missing and (b) changing `pass` to
/// `miss` if the obs value is now missing.
void assignToFloatVariable(const ufo::Variable &variable,
                           const AssignmentParameters &params,
                           const std::vector<bool> &apply,
                           const ObsFilterData &data,
                           ioda::ObsSpace &obsdb,
                           ioda::ObsDataVector<int> &qcflags) {
  ioda::ObsDataVector<float> values =
    getCurrentValues<float>(variable, obsdb, params.skipDerived);
  assignNumericValues(params, variable, apply, data, values);
  saveValues(variable, values, obsdb);
  if (variable.group() == "ObsValue" || variable.group() == "DerivedObsValue") {
    updateQCFlags(values, qcflags);
  }
}

/// Retrieve the current values of a non-numeric variable \p variable from \p obsdb (or if it
/// doesn't already exist, fill it with missing values), assign new values to elements selected by
/// the `where` clause and save the results to \p obsdb.
template <typename VariableType>
void assignToNonnumericVariable(const ufo::Variable &variable,
                                const AssignmentParameters &params,
                                const std::vector<bool> &apply,
                                const ObsFilterData &data,
                                ioda::ObsSpace &obsdb) {
  ioda::ObsDataVector<VariableType> values =
    getCurrentValues<VariableType>(variable, obsdb, params.skipDerived);
  assignNonnumericValues(params, variable, apply, data, values);
  saveValues(variable, values, obsdb);
}

/// Delegate work to an appropriate function depending on whether \p dtype is a numeric
/// or non-numeric type.
void assignToVariable(const ufo::Variable &variable,
                      ioda::ObsDtype dtype,
                      const AssignmentParameters &params,
                      const std::vector<bool> &apply,
                      const ObsFilterData &data,
                      ioda::ObsSpace &obsdb,
                      ioda::ObsDataVector<int> &qcflags) {
  switch (dtype) {
  case ioda::ObsDtype::Float:
    assignToFloatVariable(variable, params, apply, data, obsdb, qcflags);
    break;
  case ioda::ObsDtype::Integer:
    assignToIntVariable(variable, params, apply, data, obsdb);
    break;
  case ioda::ObsDtype::String:
    assignToNonnumericVariable<std::string>(variable, params, apply, data, obsdb);
    break;
  case ioda::ObsDtype::DateTime:
    assignToNonnumericVariable<util::DateTime>(variable, params, apply, data, obsdb);
    break;
  case ioda::ObsDtype::Empty:
    oops::Log::info() << "ufo::VariableAssignment::assignToVariable "
                      << "not performed on empty MPI ranks" << std::endl;
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
  if (variable.group() == "ObsValue") {
    throw eckit::BadValue("Assignment to variables from the ObsValue group is not allowed",
                          Here());
  }
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
      if (obsdb.has(variable.group(), variableWithChannel)) {
        return obsdb.dtype(variable.group(), variableWithChannel);
      }
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
  const int numOptionsSet = static_cast<int>(value_.value() != boost::none) +
                            static_cast<int>(sourceVariable.value() != boost::none) +
                            static_cast<int>(function.value() != boost::none);
  if (numOptionsSet != 1) {
    throw eckit::UserError(path.path() +
                           ": Exactly one of the 'value', 'source variable' and 'function' options "
                           "must be present");
  }
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
    if (assignment.sourceVariable.value() != boost::none) {
      allvars_ += *assignment.sourceVariable.value();
    }
    if (assignment.function.value() != boost::none) {
      allvars_ += *assignment.function.value();
    }
  }
}

void VariableAssignment::doFilter() {
  oops::Log::trace() << "VariableAssignment doFilter begin" << std::endl;

  // Select locations at which the filter will be applied
  const std::vector<bool> apply = processWhere(parameters_.where, data_, parameters_.whereOperator);

  // Assign values to successive sets of variables
  for (const AssignmentParameters &assignment : parameters_.assignments.value()) {
    const ufo::Variable variable = getVariable(assignment);
    const ioda::ObsDtype dtype = getDataType(assignment.type, variable, obsdb_);
    assignToVariable(variable, dtype, assignment, apply, data_, obsdb_, *flags_);
  }

  oops::Log::trace() << "VariableAssignment doFilter end" << std::endl;
}

void VariableAssignment::print(std::ostream & os) const {
  os << "VariableAssignment: config = " << parameters_ << std::endl;
}

}  // namespace ufo
