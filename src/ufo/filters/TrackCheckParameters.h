/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKPARAMETERS_H_
#define UFO_FILTERS_TRACKCHECKPARAMETERS_H_

#include <map>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling the operation of the track check filter.
class TrackCheckParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(TrackCheckParameters, TrackCheckUtilsParameters)

 public:
  /// Assumed temporal resolution of the observations, i.e. absolute accuracy of the reported
  /// observation times.
  oops::Parameter<util::Duration> temporalResolution{
    "temporal_resolution", util::Duration("PT1M"), this};
  /// Assumed spatial resolution of the observations (in km), i.e. absolute accuracy of the
  /// reported positions.
  ///
  /// Instantaneous speeds are estimated conservatively with the formula
  ///
  /// speed_estimate = (reported_distance - spatial_resolution) /
  ///                  (reported_time + temporal_resolution).
  oops::Parameter<float> spatialResolution{
    "spatial_resolution", 1.0f, this};
  /// Controls the size of the set of observations against which each observation is compared.
  ///
  /// Each observation O(x, t) (taken at time t and location x) is compared against the smallest
  /// set of observations O'(x', t') immediately preceding and following O(x, t) that contains
  /// \c num_distinct_buddies_per_direction earlier observations meeting the following conditions:
  /// * |t' - t| > distinct_buddy_resolution_multiplier * temporal_resolution
  /// * |x' - x| > distinct_buddy_resolution_multiplier * spatial_resolution
  /// * O' has not yet been rejected
  /// and the same number of later observations meeting the same conditions.
  oops::Parameter<int> numDistinctBuddiesPerDirection{
    "num_distinct_buddies_per_direction", 3, this};
  /// Controls the size of the set of observations against which each observation is compared.
  ///
  /// \see numDistinctBuddiesPerDirection
  oops::Parameter<int> distinctBuddyResolutionMultiplier{
    "distinct_buddy_resolution_multiplier", 3, this};

  /// Maximum allowed rate of ascent and descent (Pa/s). If not set, climb rate checks are disabled.
  oops::OptionalParameter<float> maxClimbRate{"max_climb_rate", this};
  /// Encoding of the function mapping air pressure (in Pa) to the maximum speed (in m/s)
  /// considered to be realistic.
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
  oops::Parameter<float> rejectionThreshold {
    "rejection_threshold", 0.5f, this};
  /// Name of air pressure coordinate
  oops::Parameter<std::string> pressureCoord{"pressure_coordinate",
                                             "Name of air pressure coordinate",
                                             "pressure",
                                             this};
  /// Name of air pressure group
  oops::Parameter<std::string> pressureGroup{"pressure_group",
                                             "Name of air pressure group",
                                             "ObsValue",
                                             this};

  /// Consider all observations in the track check regardless of QC flags.
  /// If this option is true then only the `where` clause has an impact on the observations used.
  /// This option should be set to `true` in order to ensure compatibility with the Met Office
  /// OPS system for aircraft processing.
  oops::Parameter<bool> ignoreExistingQCFlags
    {"ignore existing QC flags",
     "Do not consider QC flags when selecting observations to use in the track check.",
     false, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKPARAMETERS_H_
