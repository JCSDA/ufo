/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_FILTERPARAMETERSBASE_H_
#define UFO_FILTERS_FILTERPARAMETERSBASE_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"  // for FilterActionFactory
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Parameters controlling the action performed on observations flagged by a filter.
class FilterActionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FilterActionParameters, Parameters)

 public:
  /// After deserialization, holds an instance of a subclass of FilterActionParametersBase
  /// controlling the behavior of a filter action. The type of the subclass is determined
  /// by the value of the `name` key in the Configuration object from which this object
  /// is deserialized.
  oops::RequiredPolymorphicParameter<FilterActionParametersBase, FilterActionFactory>
    actionParameters{"name", this};
};


/// \brief Parameters shared by all filters derived from FilterBase.
class FilterParametersBaseWithAbstractActions : public oops::ObsFilterParametersBase {
  OOPS_ABSTRACT_PARAMETERS(FilterParametersBaseWithAbstractActions, ObsFilterParametersBase)

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

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// If set to true, the filter will be executed only after the obs operator (even if it
  /// doesn't require any variables from the GeoVaLs or HofX groups).
  oops::Parameter<bool> deferToPost{"defer to post", false, this};

  /// Return parameters specifying the actions to be performed on observations flagged by the
  /// filter.
  virtual std::vector<std::unique_ptr<FilterActionParametersBase>> actions() const = 0;

  /// \brief Parameter specifying path to yaml file containing Observation to GeoVaL name mapping
  oops::OptionalParameter<std::string> AliasFile{"observation alias file", this};

 protected:
  /// Parameters specifying a single action to be performed on observations flagged by the filter.
  ///
  /// Either `action` or `actions` may be set, but not both.
  oops::OptionalParameter<FilterActionParameters> action_{"action", this};

  /// Parameters specifying any number of actions to be performed on observations flagged by the
  /// filter.
  ///
  /// Either `action` or `actions` may be set, but not both.
  oops::OptionalParameter<std::vector<FilterActionParameters>> actions_{"actions", this};
};


/// \brief Parameters shared by all filters having a default action (typically "reject").
class FilterParametersBase : public FilterParametersBaseWithAbstractActions {
  OOPS_ABSTRACT_PARAMETERS_ENABLE_COPY_AND_MOVE(FilterParametersBase,
                                                FilterParametersBaseWithAbstractActions)

 protected:
  /// \brief Create an instance of this class configured so that if neither the `action` nor
  /// `actions` option is set, the filter performs the action \p defaultAction.
  explicit FilterParametersBase(Parameters* parent = nullptr, const char *defaultAction = "reject")
    : FilterParametersBaseWithAbstractActions(parent), defaultAction_(defaultAction)
  {}

 public:
  // Import both overloads of deserialize() from the base class. We will override one of them.
  using FilterParametersBaseWithAbstractActions::deserialize;

  /// \brief Load the values of all previously registered parameters from the (not necessarily
  /// top-level) configuration \p config.
  ///
  /// The base class implementation is overridden to throw an exception if both the `action` and
  /// `actions` keys are present in the input configuration `config`.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  /// Return parameters specifying the actions to be performed on observations flagged by the
  /// filter.
  ///
  /// This implementation returns the actions specified in the `action` or `actions` option
  /// if either is present or the default action otherwise.
  std::vector<std::unique_ptr<FilterActionParametersBase>> actions() const override;

 private:
  std::string defaultAction_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERPARAMETERSBASE_H_
