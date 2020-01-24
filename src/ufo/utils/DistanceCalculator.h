/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DISTANCECALCULATOR_H_
#define UFO_UTILS_DISTANCECALCULATOR_H_

namespace ufo
{

/// \brief Calculates distances between observations and centres of bins used during thinning.
///
/// It is assumed that the location of each point can be expressed in terms of its latitude x1,
/// longitude x2 and possibly other *nonspatial* coordinates xi (i >= 3), and that scale factors
/// si (i >= 1) are used to make all coordinates dimensionless. The distance between
/// two points is given by
///
/// distance((x1, x2, x3, x4, ...), (y1, y2, y3, y4, ...), (s1, s2, s3, s4, ...)) = finalise(
///   spatialDistanceComponent((x1, x2), (y1, y2), (s1, s2)) @
///   nonspatialDistanceComponent(x3, y3, s3) @
///   nonspatialDistanceComponent(x4, y4, s4) @ ...),
///
/// where @ is the binary operator implemented by combineDistanceComponents(). For example, to
/// use the Euclidean norm as the distance function, one would use
///
/// spatialDistanceComponent((x1, x2), (y1, y2), (s1, s2)) = (s1*(x1 - y1))**2 + (s2*(x2 - y2))**2
///                   nonspatialDistanceComponent(x, y, s) = (s*(x - y))**2
///                        combineDistanceComponents(x, y) = x + y
///                                            finalise(x) = sqrt(x).
class DistanceCalculator {
 public:
  virtual ~DistanceCalculator() {}

  virtual float spatialDistanceComponent(float obsLatitude, float obsLongitude,
                                         float latitudeBinCenter, float longitudeBinCenter,
                                         float inverseLatitudeBinWidth,
                                         float inverseLongitudeBinWidth) const = 0;

  virtual float nonspatialDistanceComponent(float obs, float binCenter,
                                            float inverseBinWidth) const = 0;

  virtual float combineDistanceComponents(float componentA, float componentB) const = 0;

  virtual float finalise(float combinedComponents) const = 0;
};

}  // namespace ufo

#endif  // UFO_UTILS_DISTANCECALCULATOR_H_
