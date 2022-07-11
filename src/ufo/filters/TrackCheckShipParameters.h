/*
 * (C) Copyright 2021 Met Office UK
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

class TrackCheckShipCoreParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TrackCheckShipCoreParameters, Parameters)

 public:
  /// Assumed temporal resolution of the observations, i.e. absolute accuracy of the reported
  /// observation times.
  oops::RequiredParameter<util::Duration> temporalResolution {
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

  /// \brief If \p earlyBreakCheck set to true, check will stop early based on the number
  /// of short-spaced, fast, and bended segments of the track
  oops::RequiredParameter<bool> earlyBreakCheck {
    "early break check", this
  };
};

/// \brief Options controlling the operation of the ship track check filter.
class TrackCheckShipParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(TrackCheckShipParameters, TrackCheckUtilsParameters)

 public:
  TrackCheckShipCoreParameters core{this};

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
  oops::Parameter<SurfaceObservationSubtype> inputCategory {
    "input category", SurfaceObservationSubtype::SHPSYN, this
  };


  /// \brief To be set to \p true if the filter's unit tests are being run
  oops::Parameter<bool> testingMode {
    "unit testing mode", false, this
  };

  /// \brief To be set to \p true if the filter's single-segment comparison test is being run.
  oops::Parameter<bool> comparisonTesting {
    "comparison test", false, this
  };

  /// Treat each record as a single observation. If this option is set to true then the records
  /// on all MPI ranks are considered together (in contrast to treating each record in isolation).
  ///
  /// The variable used to group observations into records can be set with the
  /// `obs space.obsdatain.obsgrouping.group` variable YAML option.
  oops::Parameter<bool> recordsAreSingleObs{"records_are_single_obs", false, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPPARAMETERS_H_
