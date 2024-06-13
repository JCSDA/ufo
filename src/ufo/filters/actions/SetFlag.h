/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_SETFLAG_H_
#define UFO_FILTERS_ACTIONS_SETFLAG_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

/// Observations to be skipped by the `set` or `unset` action.
enum class IgnoredObservations {
  /// No observations.
  NONE,

  /// Observations with QC flags indicating rejection.
  ///
  /// For these QC flags QCflags::isRejected() returns true.
  REJECTED,

  /// Observations with QC flags indicating irrecoverable failure.
  ///
  /// For these QC flags QCflags::isDefective() returns true.
  DEFECTIVE
};

/// Struct used to convert strings to elements of the IgnoredObservations enum.
struct IgnoredObservationsParameterTraitsHelper {
  typedef IgnoredObservations EnumType;
  static constexpr char enumTypeName[] = "IgnoredObservations";
  static constexpr util::NamedEnumerator<IgnoredObservations> namedValues[] = {
    { IgnoredObservations::NONE, "none" },
    { IgnoredObservations::REJECTED, "rejected observations" },
    { IgnoredObservations::DEFECTIVE, "defective observations" }
  };
};

}  // namespace ufo

// -----------------------------------------------------------------------------

namespace oops {

template <>
struct ParameterTraits<ufo::IgnoredObservations> :
    public EnumParameterTraits<ufo::IgnoredObservationsParameterTraitsHelper>
{};

}  // namespace oops

// -----------------------------------------------------------------------------

namespace ufo {

/// Options taken by the `set` or `unset` action.
class SetFlagParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(SetFlagParameters, FilterActionParametersBase);

 public:
  /// Name of the diagnostic flag to set or unset.
  oops::RequiredParameter<std::string> flag{"flag", this};

  /// Indicates whether the action should skip observations with certain QC flags.
  ///
  /// Allowed values:
  /// - `none` (default): Don't skip any observations.
  /// - `rejected observations`: Skip observations with QC flags indicating rejection.
  /// - `defective observations`: Skip observations with QC flags indicating irrecoverable failure
  ///   (missing observed value, rejection at the pre-processing stage or inability to compute the
  ///   model equivalent).
  oops::Parameter<IgnoredObservations> ignore{"ignore", IgnoredObservations::NONE, this};

  /// If true, this will set/unset diagnostic flags for the observation report as
  /// well as the filter variables. The observation report flag is set at a particular
  /// location if at least one filter variable has been flagged at that location.
  oops::Parameter<bool> setObservationReportFlags{"set observation report flags",
      false, this};

  /// If true, set all filter variable diagnostic flags to the value
  /// of the observation report diagnostic flag. This option can only be used if
  /// `set observation report flags` is set to true.
  oops::Parameter<bool> setFlagsToObservationReport
    {"set variable flags to observation report", false, this};
};

// -----------------------------------------------------------------------------

/// Set a diagnostic flag to \p value for all observations flagged by the filter.
template <bool value>
class SetFlag : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef SetFlagParameters Parameters_;

  explicit SetFlag(const Parameters_ &);

  void apply(const Variables & vars,
             const std::vector<std::vector<bool>> & flagged,
             ObsFilterData & data,
             int /*filterQCflag*/,
             ioda::ObsDataVector<int> & flags,
             ioda::ObsDataVector<float> & obserr) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_SETFLAG_H_
