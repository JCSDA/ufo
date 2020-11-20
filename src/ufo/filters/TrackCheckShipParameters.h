/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKSHIPPARAMETERS_H_
#define UFO_FILTERS_TRACKCHECKSHIPPARAMETERS_H_

#include <string>

#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/TrackCheckUtilsParameters.h"



namespace ufo {

/// \brief Options controlling the operation of the ship track check filter.
class TrackCheckShipParameters : public TrackCheckUtilsParameters {
 public:
  /// Assumed temporal resolution of the observations, i.e. absolute accuracy of the reported
  /// observation times.
  oops::RequiredParameter<util::Duration> temporalResolution{
    "temporal resolution", this
  };

  /// Assumed spatial resolution (km) of the observations, i.e. absolute accuracy of the
  /// reported positions.
  ///
  /// Instantaneous speeds are estimated conservatively with the formula
  ///
  /// speed_estimate = (reported_distance - spatial resolution) /
  ///                  (reported_time + temporal resolution).
  oops::RequiredParameter<double> spatialResolution {
    "spatial resolution (km)", this
  };

  /// Maximum speed (before marking as fast) in m/s
  oops::RequiredParameter<double> maxSpeed {
    "max speed (m/s)", this
  };

  /// The fraction of track observations that must be flagged in the track filter
  /// for the full track to be rejected.
  oops::RequiredParameter<float> rejectionThreshold {
    "rejection threshold", this
  };

  /// The start of an observation window where trace output should be produced. If blank,
  /// the start of the track will be treated as the start of this window.
  oops::OptionalParameter<float> debugWindowStart {
    "debug window start", this
  };

  /// The end of an observation window where trace output should be produced. If blank,
  /// the end of the track will be treated as the end of this window.
  oops::OptionalParameter<float> debugWindowEnd {
    "debug window end", this
  };

  /// The type of input source. This affects the treatment of tracks
  /// with large numbers of short segments between observations.
  oops::Parameter<int> inputCategory {
    "input category", 2, this  // 1 for buoy/other fixed input; 2 for ship
  };

  /// \brief If \p earlyBreakCheck set to true, check will stop early based on the number
  /// of short-spaced, fast, and bended segments of the track
  oops::RequiredParameter<bool> earlyBreakCheck {
    "early break check", this
  };

  /// \brief To be set to \p true if the filter's unit tests are being run
  oops::Parameter<bool> testingMode {
    "unit testing mode", false, this
  };

  /// \brief To be set to \p true if the filter's single-segment comparison test is being run.
  oops::Parameter<bool> comparisonTesting {
    "comparison test", false, this
  };
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPPARAMETERS_H_
