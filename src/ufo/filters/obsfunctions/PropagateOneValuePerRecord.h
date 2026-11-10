/*
 * (C) Copyright 2026 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_PROPAGATEONEVALUEPERRECORD_H_
#define UFO_FILTERS_OBSFUNCTIONS_PROPAGATEONEVALUEPERRECORD_H_

#include <optional>
#include <cstddef>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

enum class PropagateOneValuePerRecordPropagationDirection {
  OriginalToExtended, ExtendedToOriginal
};

struct PropagateOneValuePerRecordPropagationDirectionParameterTraitsHelper {
  typedef PropagateOneValuePerRecordPropagationDirection EnumType;
  static constexpr char enumTypeName[] = "PropagateOneValuePerRecordPropagationDirection";
  static constexpr util::NamedEnumerator<PropagateOneValuePerRecordPropagationDirection>
    namedValues[] = {
    { PropagateOneValuePerRecordPropagationDirection::OriginalToExtended, "from original" },
    { PropagateOneValuePerRecordPropagationDirection::ExtendedToOriginal, "from extended" }
  };
};
}  // namespace ufo

namespace oops {
template <>
struct ParameterTraits<ufo::PropagateOneValuePerRecordPropagationDirection> :
    public EnumParameterTraits<
      ufo::PropagateOneValuePerRecordPropagationDirectionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {


/// \brief Options controlling the PropagateOneValuePerRecord ObsFunction
template <typename PropagateVariableType>
class PropagateOneValuePerRecordParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PropagateOneValuePerRecordParameters, Parameters)

 public:
  /// If set to false, the selected value in each record is propagated to all
  /// locations in that record, but not to the companion record.
  oops::Parameter<bool> propagateToCompanionAsWell{
    "propagate to companion as well", false, this};
  /// If set to "from original", the value selected by `in-record location of value to
  /// propagate` in each original record is propagated to all locations in both the original
  /// and companion extended record. If set to "from extended", the selected value in each
  /// extended record is propagated to all locations in both the extended and companion
  /// original record. Only used if `propagate to companion as well` is true.
  oops::OptionalParameter<PropagateOneValuePerRecordPropagationDirection> propagationDirection{
     "propagation direction", this};

  /// Name of the variable whose selected value in each record is propagated to all locations in
  /// that record and, when requested, its companion.
  oops::RequiredParameter<Variable> variableToPropagate{"variable to propagate", this};

  /// If true, derived variables won't be prioritised when retrieving the variable to propagate.
  oops::Parameter<bool> skipDerived{"skip derived", false, this};

  /// In-record location index which selects the value to be propagated. Can be
  /// positive or negative, the latter means counting from the end of the record (e.g.
  /// -1 means the last location in the record, -2 means the second to last, and so on.)
  oops::RequiredParameter<int> inRecordLocationOfValueToPropagate{
     "in-record location of value to propagate", this};
};

// -----------------------------------------------------------------------------

/// \brief Propagate variables from the original ObsSpace to the extended ObsSpace or vice versa.
///
/// This ObsFunction propagates the value selected by `in-record location of value to
/// propagate` from each record to all locations in that record. When requested, it also
/// propagates the selected value to all locations in the corresponding companion record
/// in the extended (or original) ObsSpace.
template <typename PropagateVariableType>
class PropagateOneValuePerRecord : public ObsFunctionBase<PropagateVariableType> {
 public:
  explicit PropagateOneValuePerRecord(const eckit::LocalConfiguration &);
  ~PropagateOneValuePerRecord();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<PropagateVariableType> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  PropagateOneValuePerRecordParameters<PropagateVariableType> options_;
  ufo::Variables invars_;
  std::optional<std::size_t> findLocIndex(const int value,
     const std::size_t locsSize) const;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_PROPAGATEONEVALUEPERRECORD_H_
