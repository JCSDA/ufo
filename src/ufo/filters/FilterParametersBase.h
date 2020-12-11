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
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

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

  /// Operation performed on observations meeting the conditions specified in the `where` clause.
  oops::Parameter<eckit::LocalConfiguration> action{"action", eckit::LocalConfiguration(), this};

  /// If set to true, the filter will be executed only after the obs operator (even if it
  /// doesn't require any variables from the GeoVaLs or HofX groups).
  oops::Parameter<bool> deferToPost{"defer to post", false, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERPARAMETERSBASE_H_
