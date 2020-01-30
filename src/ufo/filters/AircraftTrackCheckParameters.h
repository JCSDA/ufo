/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_AIRCRAFTTRACKCHECKPARAMETERS_H_
#define UFO_FILTERS_AIRCRAFTTRACKCHECKPARAMETERS_H_

#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraitsMap.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// \brief Options controlling the operation of the aircraft track check filter.
class AircraftTrackCheckParameters : public oops::Parameters {
 public:
  /// Assumed temporal resolution of the observations, i.e. absolute accuracy of the reported
  /// observation times.
  oops::Parameter<util::Duration> temporalResolution{
    "temporal_resolution", util::Duration("PT1M"), this};
  /// Assumed spatial resolution of the observations (in km), i.e. absolute accuracy of the
  /// reported aircraft positions.
  ///
  /// Instantaneous speeds are estimated conservatively with the formula
  ///
  /// speed_estimate = (reported_distance - spatial_resolution) /
  ///                  (reported_time + temporal_resolution).
  oops::Parameter<float> spatialResolution{
    "spatial_resolution", 1.0f, this};

  /// Controls the size of the set of observations against which each observation is compared.
  ///
  /// Each observation O(x, t) (taken at time t and location x) is compared against the set of all
  /// observations O'(x', t') immediately preceding and following it that satisfy at least one of
  /// the following two conditions:
  /// * |t - t'| <= core_temporal_neighborhood_radius
  /// * |x - x'| <= core_spatial_neighborhood_radius
  /// and against half_num_noncore_neighbors non-rejected observations immediately following that
  /// set and the same number of non-rejected observations immediately preceding that set.
  oops::Parameter<util::Duration> coreTemporalNeighborhoodRadius{
    "core_temporal_neighborhood_radius", util::Duration("PT2M"), this};
  /// Controls the size of the set of observations against which each observation is compared.
  ///
  /// \see coreTemporalNeighborhoodRadius
  oops::Parameter<float> coreSpatialNeighborhoodRadius{
    "core_spatial_neighborhood_radius", 2.0f, this};
  /// Controls the size of the set of observations against which each observation is compared.
  ///
  /// \see coreTemporalNeighborhoodRadius
  oops::Parameter<int> halfNumNoncoreNeighbors{
    "half_num_noncore_neighbors", 3, this};

  /// Maximum allowed rate of ascent and descent (Pa/s).
  oops::Parameter<float> maxClimbRate{"max_climb_rate", 300.0f, this};
  /// Encoding of the function mapping air pressure (in Pa) to the maximum allowed aircraft speed
  /// (in m/s).
  ///
  /// The function is taken to be a linear interpolation of a series of (pressure, speed) points.
  /// The pressures and speeds at these points should be specified as keys and values of a
  /// JSON-style map. Owing to a bug in the eckit YAML parser, the keys must be enclosed in quotes.
  /// For example,
  ///
  ///   max_speed_interpolation_points: { "0": 900, "100000": 100 }
  ///
  /// encodes a linear function equal to 900 m/s at 0 Pa and 100 m/s at 100000 Pa.
  oops::Parameter<std::map<float, float>> maxSpeedInterpolationPoints{
    "max_speed_interpolation_points", std::map<float, float>{{0.0f, 1000.0f}}, this};
  /// Maximum fraction of climb rate or speed estimates obtained by comparison with other
  /// observations that are allowed to fall outside the allowed ranges before an observation is
  /// rejected.
  oops::Parameter<float> rejectionThreshold{
    "rejection_threshold", 0.5f, this};

  /// Variable storing integer-valued flight IDs. Observations taken during each flight are
  /// checked separately.
  oops::Parameter<Variable> flightIdVariable{
    "flight_id_variable", Variable("flight_id@MetaData"), this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_AIRCRAFTTRACKCHECKPARAMETERS_H_
