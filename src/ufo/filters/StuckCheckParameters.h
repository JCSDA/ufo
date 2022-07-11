/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_STUCKCHECKPARAMETERS_H_
#define UFO_FILTERS_STUCKCHECKPARAMETERS_H_

#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"

namespace ufo {
class StuckCheckCoreParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StuckCheckCoreParameters, Parameters)

 public:
  /// The maximum number of observations in a row that can have the same observed value
  /// before the observations may be flagged by this filter.
  oops::OptionalParameter<size_t> numberStuckTolerance {
    "number stuck tolerance", this
  };

  /// The maximum duration in which an observation value is "stuck" before
  /// the observations may be flagged (unless all of the observations have that one value,
  /// in which case the observations could be flagged anyway)
  oops::OptionalParameter<util::Duration> timeStuckTolerance {
    "time stuck tolerance", this
  };

  /// The maximum percentage of observations in a row, out of the total number of observations
  /// in the record, that can have the same observed value before the observations may be
  /// flagged by this filter.
  oops::OptionalParameter<float> percentageStuckTolerance {
    "percentage stuck tolerance", this
  };
};

/// \brief Options controlling the operation of stuck check filter.
class StuckCheckParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(StuckCheckParameters, TrackCheckUtilsParameters)
 public:
  StuckCheckCoreParameters core{this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_STUCKCHECKPARAMETERS_H_
