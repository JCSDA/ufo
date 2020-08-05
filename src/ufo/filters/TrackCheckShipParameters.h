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

#include "ufo/filters/TrackCheckUtilsParameters.h"



namespace ufo {

/// \brief Options controlling the operation of the ship track check filter. TODO: edit
class TrackCheckShipParameters : public TrackCheckUtilsParameters {
 public:
  /// Assumed temporal resolution of the observations, i.e. absolute accuracy of the reported
  /// observation times.
  oops::Parameter<util::Duration> temporalResolution{
    "temporal resolution", util::Duration("PT10M"), this
  };
  /// Assumed spatial resolution of the observations (in km), i.e. absolute accuracy of the
  /// reported positions.
  ///
  /// Instantaneous speeds are estimated conservatively with the formula
  ///
  /// speed_estimate = (reported_distance - spatial resolution) /
  ///                  (reported_time + temporal resolution).
  oops::Parameter<double> spatialResolution {
    "spatial resolution", 1.0, this
  };

  /// Maximum speed (before marking as fast) in km/s
  /// \attention: please let me know what a better default value for ships might be.
  oops::Parameter<double> maxSpeed {
    "max speed", 1.0, this
  };

  oops::OptionalParameter<float> debugWindowStart {
    "debug window start", this
  };

  oops::OptionalParameter<float> debugWindowEnd {
    "debug window end", this
  };

  enum class InputCategory : int
  {
    FIXED = 1,
    SHIP = 2,
    SHIP_AUTO = 3  // GPS unavailable; ship assumed
  };

  oops::Parameter<int> inputCategory {
    "input category", 2, this  // 1 for buoy/other fixed input; 2 for ship; 3 for ship_auto
  };

  /// \brief If \p deferredCheckSimultaneous set to true, check of simultaneous
  /// observations will be deferred until the end, at which point flagged observations
  /// may be reinstated.
  ///
  /// \todo implementation in progress
  oops::Parameter<bool> deferredCheckSimultaneous {
    "deferred check simultaneous", false, this
  };

  /// \brief If \p earlyBreakCheck set to true, check will stop early based on the number
  /// of short-spaced, fast, and bended segments of the track
  ///
  /// \todo implementation in progress
  oops::Parameter<bool> earlyBreakCheck {
    "early break check", true, this
  };

  /// To be set to \p true if the filter's tests are being run
  oops::Parameter<bool> testingMode {
    "testing mode", false, this
  };
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPPARAMETERS_H_
