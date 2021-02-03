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
  /// is deserialized. Omitting this key is equivalent to setting it to `reject`.
  oops::PolymorphicParameter<FilterActionParametersBase, FilterActionFactory>
    actionParameters{"name", "reject", this};
};

/// \brief Parameters shared by all filters derived from FilterBase.
class FilterParametersBase : public oops::ObsFilterParametersBase {
  OOPS_ABSTRACT_PARAMETERS(FilterParametersBase, ObsFilterParametersBase)

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
  oops::Parameter<eckit::LocalConfiguration> where{"where", eckit::LocalConfiguration(), this};

  /// Parameters controlling the action performed on observations flagged by the filter.
  oops::Parameter<FilterActionParameters> action{"action", {}, this};

  /// If set to true, the filter will be executed only after the obs operator (even if it
  /// doesn't require any variables from the GeoVaLs or HofX groups).
  oops::Parameter<bool> deferToPost{"defer to post", false, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERPARAMETERSBASE_H_
