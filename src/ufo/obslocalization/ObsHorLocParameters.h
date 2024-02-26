/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_

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
  GEODESIC, CARTESIAN
};

enum class SearchMethod {
  BRUTEFORCE, KDTREE
};

struct DistanceTypeParameterTraitsHelper {
  typedef DistanceType EnumType;
  static constexpr char enumTypeName[] = "DistanceType";
  static constexpr util::NamedEnumerator<DistanceType> namedValues[] = {
    { DistanceType::GEODESIC, "geodesic" },
    { DistanceType::CARTESIAN, "cartesian" }
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

  /// Distance calculation type: geodesic (on sphere) or cartesian (euclidian)
  /// Default: geodesic
  oops::Parameter<DistanceType> distanceType{"distance type", DistanceType::GEODESIC, this};

  /// returns distance between points \p p1 and \p p2, depending on the
  /// distance calculation type distanceType
  double distance(const eckit::geometry::Point3 & p1, const eckit::geometry::Point3 & p2) const {
    if (distanceType == DistanceType::GEODESIC) {
      eckit::geometry::Point2 q1(p1[0], p1[1]);
      eckit::geometry::Point2 q2(p1[0], p1[1]);
      return eckit::geometry::Sphere::distance(radius_earth, q1, q2);
    } else {
      ASSERT(distanceType == DistanceType::CARTESIAN);
      return p1.distance(p2);
    }
  }

  // Earth radius in m
  static constexpr double radius_earth = 6.371e6;
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCPARAMETERS_H_
