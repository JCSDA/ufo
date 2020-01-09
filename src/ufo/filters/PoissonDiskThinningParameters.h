/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_POISSONDISKTHINNINGPARAMETERS_H_
#define UFO_FILTERS_POISSONDISKTHINNINGPARAMETERS_H_

#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

enum class ExclusionVolumeShape {
  CYLINDER, ELLIPSOID
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::ExclusionVolumeShape> {
  static boost::optional<ufo::ExclusionVolumeShape> get(const eckit::Configuration &config,
                                                        const std::string& name) {
    std::string value;
    if (config.get(name, value)) {
      if (value == "cylinder")
        return ufo::ExclusionVolumeShape::CYLINDER;
      if (value == "ellipsoid")
        return ufo::ExclusionVolumeShape::ELLIPSOID;
      throw eckit::BadParameter("Bad conversion from std::string '" + value +
                                "' to ExclusionVolumeShape", Here());
    } else {
      return boost::none;
    }
  }
};

}  // namespace oops

namespace ufo {

/// \brief Options controlling the operation of the PoissonDiskThinning filter.
///
/// \note The descriptions of several options refer to the _exclusion volume_, which is a domain
/// surrounding the location of each observation. If an observation is retained, then no other
/// observations lying in the interior of its exclusion volume may be retained at the same time.
class PoissonDiskThinningParameters : public oops::Parameters {
 public:
  // Exclusion volume

  /// Size of the exclusion volume in the horizontal direction (in km).
  ///
  /// If this size needs to vary from observation to observation, this parameter should be set
  /// to the name of a variable storing the exclusion volume size for each observation. Otherwise
  /// the parameter can be set to floating-point number. If the parameter is not set, horizontal
  /// position is ignored during thinning.
  oops::OptionalParameter<std::string> minHorizontalSpacing{"min_horizontal_spacing", this};

  /// Size of the exclusion volume in the vertical direction (in Pa).
  ///
  /// Like min_horizontal_spacing, this can be either a variable name or a floating-point number.
  oops::OptionalParameter<std::string> minVerticalSpacing{"min_vertical_spacing", this};

  /// Size of the exclusion volume in the temporal direction.
  ///
  /// This must be a duration (not a variable name).
  // TODO(wsmigaj): allow this option to be set to a variable name. See comment in
  // PoissonDiskThinning::getObsData() for explanation why this isn't currently supported.
  oops::OptionalParameter<util::Duration> minTimeSpacing{"min_time_spacing", this};

  /// Shape of the exclusion volume surrounding each observation.
  ///
  /// Allowed values:
  /// - \c cylinder: the exclusion volume of an observation taken at latitude lat, longitude lon,
  ///   pressure p and time t is the set of all locations (lat', lon', p', t') for which all of
  ///   the following conditions are met:
  ///   * the geodesic distance between (lat, lon) and (lat', lon') is smaller than
  ///     min_horizontal_spacing
  ///   * |p - p'| < min_vertical_spacing
  ///   * |t - t'| < min_time_spacing.
  /// - \c ellipsoid: the exclusion volume of an observation taken at latitude lat, longitude lon,
  ///   pressure p and time t is the set of all locations (lat', lon', p', t') for which
  ///   the following condition is met:
  ///   * geodesic_distance((lat, lon), (lat', lon'))^2 / min_horizontal_spacing^2 +
  ///     (p - p')^2 / min_vertical_spacing^2 + (t - t')^2 / min_time_spacing^2 < 1.
  oops::Parameter<ExclusionVolumeShape> exclusionVolumeShape{"exclusion_volume_shape",
                                                             ExclusionVolumeShape::CYLINDER, this};

  // Observation categories

  /// Variable storing integer-valued IDs associated with observations. Observations belonging
  /// to different categories are thinned separately.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};

  // Selection of observations to retain

  /// Variable storing observation priorities. An observation will not be retained if it lies
  /// within the exclusion volume of an observation with a higher priority.
  ///
  /// \note The implementation assumes that the exclusion volumes of all observations with the same
  /// priority have the same size, and the exclusion volume size increases or stays the same with
  /// decreasing priority.
  ///
  /// If this parameter is not set, all observations are assumed to have equal priority.
  oops::OptionalParameter<Variable> priorityVariable{"priority_variable", this};

  /// If true, observations will be randomly shuffled before being inspected as candidates
  /// for retaining.
  ///
  /// \note It is recommended to leave shuffling enabled in production code, since the performance
  /// of the spatial point index (kd-tree) used in the filter's implementation may be degraded if
  /// observation locations are ordered largely monotonically (and random shuffling essentially
  /// prevents that from happening).
  oops::Parameter<bool> shuffle{"shuffle", true, this};

  /// Seed with which to initialize the random number generator used to shuffle the observations
  /// if \p shuffle is set to true.
  ///
  /// If omitted, a seed will be generated based on the current (calendar) time.
  oops::OptionalParameter<int> randomSeed{"random_seed", this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_POISSONDISKTHINNINGPARAMETERS_H_
