/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_FILTERPARAMETERSBASE_H_
#define UFO_FILTERS_FILTERPARAMETERSBASE_H_

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"  // for FilterActionFactory
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Parameters controlling the action performed on observations flagged by a filter.
class FilterActionParameters : public oops::Parameters {
  // We call this macro instead of "plain" OOPS_CONCRETE_PARAMETERS() because we want to customize
  // the constructor definition (adding an optional parameter `defaultAction`).
  OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(FilterActionParameters, Parameters)

 public:
  /// Create an instance of this class that by default (i.e. unless overridden by settings loaded
  /// from a Configuration object during deserialization) represents the action \p defaultAction.
  explicit FilterActionParameters(Parameters *parent = nullptr,
                                  const char *defaultAction = "reject")
    : Parameters(parent), actionParameters("name", defaultAction, this)
  {}

  /// After deserialization, holds an instance of a subclass of FilterActionParametersBase
  /// controlling the behavior of a filter action. The type of the subclass is determined
  /// by the value of the `name` key in the Configuration object from which this object
  /// is deserialized.
  oops::PolymorphicParameter<FilterActionParametersBase, FilterActionFactory> actionParameters;
};

/// \brief Parameters shared by all filters derived from FilterBase.
class FilterParametersBaseWithAbstractAction : public oops::ObsFilterParametersBase {
  OOPS_ABSTRACT_PARAMETERS(FilterParametersBaseWithAbstractAction, ObsFilterParametersBase)

 public:
  /// Variables (and channels) to which the filter should be applied. These will be the only
  /// variables whose QC flags or error estimates may be modified. If not specified, defaults to
  /// all simulated variables in the ObsSpace.
  oops::OptionalParameter<std::vector<Variable>> filterVariables{
    "filter variables", this};

  /// Conditions used to select observations to which the filter should be applied.
  /// If not specified, the filter will be applied to all observations.
  ///
  /// \note The QC flags of observations to which a filter is not applied won't normally be
  /// modified (i.e. those observations will automatically "pass", unless of course they have
  /// already been rejected by another filter). The only exception is the Domain Check filter,
  /// which does the exact opposite: it rejects all observations that haven't been selected by the
  /// `where` statement.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// If set to true, the filter will be executed only after the obs operator (even if it
  /// doesn't require any variables from the GeoVaLs or HofX groups).
  oops::Parameter<bool> deferToPost{"defer to post", false, this};

  /// Return parameters defining the action performed on observations flagged by the filter.
  virtual const FilterActionParametersBase &action() const = 0;
};

/// \brief Parameters shared by all filters having a default action (typically "reject").
class FilterParametersBase : public FilterParametersBaseWithAbstractAction {
  OOPS_ABSTRACT_PARAMETERS_ENABLE_COPY_AND_MOVE(FilterParametersBase,
                                                FilterParametersBaseWithAbstractAction)

 protected:
  /// \brief Create an instance of this class configured so that if the filter action is not
  /// specified, it defaults to \p defaultAction.
  explicit FilterParametersBase(Parameters* parent = nullptr, const char *defaultAction = "reject")
    : FilterParametersBaseWithAbstractAction(parent),
      action_("action", FilterActionParameters(nullptr, defaultAction), this)
  {}

 public:
  const FilterActionParametersBase &action() const override {
    return action_.value().actionParameters.value();
  }

 private:
  /// Parameters controlling the action performed on observations flagged by the filter.
  oops::Parameter<FilterActionParameters> action_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERPARAMETERSBASE_H_
