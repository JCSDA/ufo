/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_TEMPORALTHINNINGPARAMETERS_H_
#define UFO_FILTERS_TEMPORALTHINNINGPARAMETERS_H_

#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// \brief Options controlling the operation of the TemporalThinning filter.
class TemporalThinningParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TemporalThinningParameters, Parameters)

 public:
  /// Minimum spacing between two successive retained observations.
  oops::Parameter<util::Duration> minSpacing{"min_spacing", util::Duration("PT1H"), this};

  /// Only relevant if \c priority_variable is set.
  ///
  /// If \c tolerance is nonzero, then whenever an observation O lying at least \c min_spacing
  /// from the previous retained observation O' is found, the filter will inspect all observations
  /// lying no more than \c tolerance further from O' and retain the one with the highest priority.
  oops::Parameter<util::Duration> tolerance{"tolerance", util::Duration("PT0H"), this};

  /// If not set, the thinning filter will consider observations as candidates for retaining
  /// in chronological order.
  ///
  /// If set, the filter will start from the observation taken as close as possible to \c seed_time,
  /// then consider all successive observations in chronological order, and finally all preceding
  /// observations in reverse chronological order.
  oops::OptionalParameter<util::DateTime> seedTime{"seed_time", this};

  /// A string- or integer-valued variable. Observations with different values of that variable will
  /// be thinned separately.
  ///
  /// If not set and observations were grouped into records when the observation space was
  /// constructed, observations from each record will be thinned separately. If not set and
  /// observations were not grouped into records, all observations will be thinned together.
  ///
  /// Note: the variable used to group observations into records can be set with the
  /// \c obs space.obsdatain.obsgrouping.group variable YAML option.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};

  /// Variable storing observation priorities. Used together with \c tolerance; see the
  /// documentation of that parameter for more information.
  oops::OptionalParameter<Variable> priorityVariable{"priority_variable", this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_TEMPORALTHINNINGPARAMETERS_H_
