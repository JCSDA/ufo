/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_

#include <cmath>
#include <string>
#include <utility>

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/UnitSphere.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {

enum class DistanceType {
  GEODESIC, CARTESIAN, PERIODIC
};

enum class SearchMethod {
  BRUTEFORCE, KDTREE
};

struct DistanceTypeParameterTraitsHelper {
  typedef DistanceType EnumType;
  static constexpr char enumTypeName[] = "DistanceType";
  static constexpr util::NamedEnumerator<DistanceType> namedValues[] = {
    { DistanceType::GEODESIC, "geodesic" },
    { DistanceType::CARTESIAN, "cartesian" },
    { DistanceType::PERIODIC, "periodic" }
  };
};

struct SearchMethodParameterTraitsHelper {
  typedef SearchMethod EnumType;
  static constexpr char enumTypeName[] = "SearchMethod";
  static constexpr util::NamedEnumerator<SearchMethod> namedValues[] = {
    { SearchMethod::BRUTEFORCE, "brute_force" },
    { SearchMethod::KDTREE, "kd_tree" }
  };
};

}  // namespace ufo

namespace oops {

/// Extraction of DistanceType parameters from config
template <>
struct ParameterTraits<ufo::DistanceType> :
    public EnumParameterTraits<ufo::DistanceTypeParameterTraitsHelper>
{};


/// Extraction of SearchMethod parameters from config
template <>
struct ParameterTraits<ufo::SearchMethod> :
    public EnumParameterTraits<ufo::SearchMethodParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Options controlling local observations subsetting
class ObsHorLocParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsHorLocParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> method{"localization method", this};

  /// Localization lengthscale (find all obs within the distance from reference point)
  oops::Parameter<double> lengthscale{"lengthscale", 0.0, this};

  /// Method for searching for nearest points: brute force or KD-tree
  oops::Parameter<SearchMethod> searchMethod{"search method", SearchMethod::KDTREE, this};

  /// Maximum number of obs
  oops::OptionalParameter<int> maxnobs{"max nobs", this};

  /// Distance calculation type: geodesic (on sphere), cartesian (euclidian) or periodic
  /// (planar distance on a doubly periodic domain, see longitude and latitude periods below)
  /// Default: geodesic
  oops::Parameter<DistanceType> distanceType{"distance type", DistanceType::GEODESIC, this};

  /// Periods (in degrees) of a periodic domain in longitude and latitude, used with the periodic
  /// distance type. The domain is treated as planar: the distance between two points is the earth
  /// radius times the euclidian norm of their longitude and latitude differences (in radians),
  /// each difference being taken modulo its period (nearest periodic image). The longitude and
  /// latitude of a cartesian domain should be defined accordingly, as its x and y coordinates
  /// divided by the earth radius. Without a period, there is no wrap-around in that direction.
  oops::OptionalParameter<double> longitudePeriod{"longitude period", this};
  oops::OptionalParameter<double> latitudePeriod{"latitude period", this};

  /// returns distance between points \p p1 and \p p2, depending on the
  /// distance calculation type distanceType
  double distance(const eckit::geometry::Point3 & p1, const eckit::geometry::Point3 & p2) const {
    if (distanceType == DistanceType::GEODESIC) {
      eckit::geometry::Point2 q1(p1[0], p1[1]);
      eckit::geometry::Point2 q2(p2[0], p2[1]);
      return eckit::geometry::Sphere::distance(radius_earth, q1, q2);
    } else if (distanceType == DistanceType::PERIODIC) {
      const double dlon = periodicDifference(p1[0]-p2[0], longitudePeriod.value());
      const double dlat = periodicDifference(p1[1]-p2[1], latitudePeriod.value());
      return radius_earth*M_PI/180.0*std::sqrt(dlon*dlon+dlat*dlat);
    } else {
      ASSERT(distanceType == DistanceType::CARTESIAN);
      return p1.distance(p2);
    }
  }

  /// returns the difference \p diff of two coordinates, wrapped to the nearest periodic image
  /// if a period is given
  static double periodicDifference(const double diff, const boost::optional<double> & period) {
    if (period == boost::none) return diff;
    ASSERT(*period > 0.0);
    return diff - *period*std::round(diff / *period);
  }

  // Earth radius in m
  static constexpr double radius_earth = 6.371e6;
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_
