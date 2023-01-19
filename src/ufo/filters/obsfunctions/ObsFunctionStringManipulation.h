/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSTRINGMANIPULATION_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSTRINGMANIPULATION_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"


namespace ufo {

enum class StringManipOption {
  STRINGCUT
};

struct StringManipOptionParameterTraitsHelper {
  typedef StringManipOption EnumType;
  static constexpr char enumTypeName[] = "StringManipOption";
  static constexpr util::NamedEnumerator<StringManipOption> namedValues[] = {
    { StringManipOption::STRINGCUT, "stringcut" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::StringManipOption> :
    public EnumParameterTraits<ufo::StringManipOptionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsFunctionStringManipulation ObsFunction
class StringManipulationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StringManipulationParameters, Parameters)

 public:
  /// This defines the string manipulation to be carried out
  oops::RequiredParameter<StringManipOption> stringManipOption{"string operation", this};

  /// Input variable containing the string to be manipulated
  oops::RequiredParameter<std::vector<Variable>> variable{"variable", this};

  /// minimum index to start the string cut from. Starting with 0 being the far left of a string.
  oops::OptionalParameter<int> startIndex{"startIndex", this};

  /// Length of the string to be cut
  oops::OptionalParameter<int> cutLength{"cutLength", this};
};

// -----------------------------------------------------------------------------

/// \brief Tools for string manipulation
///
/// Output a string cut from an index to a set length. For example a
/// site name with the form AAAA-BBBB could be cut to AAAA by setting
/// startIndex to 0 and cutLength to 4.
///
/// Example 1
///   - filter: Variable Assignment
///     assignments:
///     - name: MetaData/stationName
///       type: string
///       function:
///         name: StringObsFunction/StringManipulation
///         options:
///         string operation: stringcut
///         variable: [MetaData/stationACCombined]
///         startIndex: 0
///         cutLength: 4

// -----------------------------------------------------------------------------

class StringManipulation : public ObsFunctionBase<std::string> {
 public:
  explicit StringManipulation(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<std::string> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  StringManipulationParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSTRINGMANIPULATION_H_
