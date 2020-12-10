/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_POISSONDISKTHINNINGPARAMETERS_H_
#define UFO_FILTERS_POISSONDISKTHINNINGPARAMETERS_H_

#include <string>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraitsScalarOrMap.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

enum class ExclusionVolumeShape {
  CYLINDER, ELLIPSOID
};

struct ExclusionVolumeShapeParameterTraitsHelper {
  typedef ExclusionVolumeShape EnumType;
  static constexpr char enumTypeName[] = "ExclusionVolumeShape";
  static constexpr util::NamedEnumerator<ExclusionVolumeShape> namedValues[] = {
    { ExclusionVolumeShape::CYLINDER, "cylinder" },
    { ExclusionVolumeShape::ELLIPSOID, "ellipsoid" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::ExclusionVolumeShape> :
    public EnumParameterTraits<ufo::ExclusionVolumeShapeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Options controlling the operation of the PoissonDiskThinning filter.
///
/// \note The descriptions of several options refer to the _exclusion volume_, which is a domain
/// surrounding the location of each observation. If an observation is retained, then no other
/// observations lying in the interior of its exclusion volume may be retained at the same time.
class PoissonDiskThinningParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PoissonDiskThinningParameters, Parameters)

 public:
  typedef int Priority;

  // Exclusion volume

  /// Size of the exclusion volume in the horizontal direction (in km).
  ///
  /// If the priority_variable parameter is not set and hence all observations have the same
  /// priority, this parameter must be a floating-point constant. Otherwise it may also be a map
  /// assigning an exclusion volume size to each observation priority. Exclusion volumes of
  /// lower-priority observations must be at least as large as those of higher-priority ones.
  /// If this parameter is not set, horizontal position is ignored during thinning.
  ///
  /// \note Owing to a bug in the eckit YAML parser, maps need to be written in the JSON style,
  /// with keys quoted. Example:
  ///
  ///   min_horizontal_spacing: {"1": 123, "2": 321}
  ///
  /// This will not work:
  ///
  ///   min_horizontal_spacing: {1: 123, 2: 321}
  ///
  /// and neither will this:
  ///
  ///   min_horizontal_spacing:
  ///     1: 123
  ///     2: 321
  ///
  /// or this:
  ///
  ///   min_horizontal_spacing:
  ///     "1": 123
  ///     "2": 321
  oops::OptionalParameter<util::ScalarOrMap<Priority, float>> minHorizontalSpacing{
    "min_horizontal_spacing", this};

  /// Size of the exclusion volume in the vertical direction (in Pa).
  ///
  /// Like min_horizontal_spacing, this can be either a constant or a map.
  oops::OptionalParameter<util::ScalarOrMap<Priority, float>> minVerticalSpacing{
    "min_vertical_spacing", this};

  /// Size of the exclusion volume in the temporal direction.
  ///
  /// Like min_horizontal_spacing, this can be either a constant or a map.
  oops::OptionalParameter<util::ScalarOrMap<Priority, util::Duration>> minTimeSpacing{
    "min_time_spacing", this};

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

  /// A string-valued or integer-valued variable. Observations with different values of that
  /// variable will be thinned separately.
  ///
  /// If not set and observations were grouped into records when the observation space was
  /// constructed, observations from each record will be thinned separately. If not set and
  /// observations were not grouped into records, all observations will be thinned together.
  ///
  /// Note: the variable used to group observations into records can be set with the
  /// \c obs space.obsdatain.obsgrouping.group variable YAML option.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};

  // Selection of observations to retain

  /// Variable storing observation priorities. An observation will not be retained if it lies
  /// within the exclusion volume of an observation with a higher priority.
  ///
  /// As noted in the documentation of min_horizontal_spacing, the exclusion volume size must be a
  /// (weakly) monotonically decreasing function of observation priority, i.e. the exclusion volumes
  /// of all observations with the same priority must have the same size, and the exclusion volumes
  /// of lower-priority observations must be at least as large as those of higher-priority ones.
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
