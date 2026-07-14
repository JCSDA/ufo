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
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
class StuckCheckCoreParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StuckCheckCoreParameters, Parameters)

 public:
  /// The maximum number of observations in a row that can have the same observed value
  /// before the observations may be flagged by this filter.
  oops::OptionalParameter<size_t> numberStuckTolerance {
    "number stuck tolerance", this
  };

  /// A non-negative integer variable containing the number stuck tolerance. All
  /// non-missing values in a record must be equal. If all values are missing
  /// then the stuck check will not be performed for that record.
  oops::OptionalParameter<Variable> numberStuckToleranceVariable {
    "number stuck tolerance variable", this
  };

  /// The maximum duration in which an observation value is "stuck" before
  /// the observations may be flagged (unless all of the observations have that one value,
  /// in which case the observations could be flagged anyway)
  oops::OptionalParameter<util::Duration> timeStuckTolerance {
    "time stuck tolerance", this
  };

  /// The maximum percentage tolerance for stuck observations - number stuck tolerance is
  /// derived from the given percentage. The calculation basis depends on whether "percentage stuck
  /// tolerance" is true or false. If false (default), the number stuck tolerance is
  ///   round(percentageStuckTolerance / 100 * num_observations_in_record)
  /// If true, the number stuck tolerance is
  ///   round(percentageStuckTolerance / 100 * num_intervals_between_observations_in_record)
  /// where num_intervals_between_observations_in_record is one less than num_observations_in_record
  oops::OptionalParameter<float> percentageStuckTolerance {
    "percentage stuck tolerance", this
  };

  /// Where the percentage stuck tolerance is set, a minimum number of stuck observations
  /// that must be exceeded before observations may be flagged. A value of 1 means that
  /// at least 2 identical observations are required to form a streak. A value of 2 means
  /// that at least 3 identical observations are required, and so on.
  oops::Parameter<size_t> minimumAllowedStuck {
    "minimum allowed stuck", 1, this
  };

  /// When true, percentageStuckTolerance is computed as a percentage of the number of
  /// consecutive identical values (intervals). When false (default), it is computed as a
  /// percentage of the total observations in the record. See "percentage stuck tolerance" for more
  /// details.
  oops::Parameter<bool> percentageStuckToleranceBasedOnIntervals {
    "percentage stuck tolerance based on intervals", false, this
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
