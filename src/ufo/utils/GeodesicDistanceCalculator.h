/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_GEODESICDISTANCECALCULATOR_H_
#define UFO_UTILS_GEODESICDISTANCECALCULATOR_H_

#include <cmath>

#include "ufo/utils/Constants.h"
#include "ufo/utils/DistanceCalculator.h"

namespace ufo
{

/// A DistanceCalculator implementing the distance function that maps each pair of points to
/// the shortest path on a spherical Earth connecting these points. (All coordinates except
/// the latitude and longitude are ignored.)
class GeodesicDistanceCalculator : public DistanceCalculator {
 public:
  float spatialDistanceComponent(float obsLatitude, float obsLongitude,
                                 float latitudeBinCenter, float longitudeBinCenter,
                                 float /*inverseLatitudeBinWidth*/,
                                 float /*inverseLongitudeBinWidth*/) const override {
    const float deg2rad = static_cast<float>(Constants::deg2rad);
    const float re = static_cast<float>(Constants::mean_earth_rad);  // km

    float q1 = std::cos((obsLongitude - longitudeBinCenter) * deg2rad);
    float q2 = std::cos((obsLatitude - latitudeBinCenter) * deg2rad);
    float q3 = std::cos((obsLatitude + latitudeBinCenter) * deg2rad);

    float geodesicDistance = (re * std::acos(0.5f*((1.0f+q1)*q2 - (1.0f-q1)*q3)) + 1.0f);
    return geodesicDistance;
  }

  float nonspatialDistanceComponent(float /*obs*/, float /*binCenter*/,
                                    float /*inverseBinWidth*/) const override {
    return 0;
  }

  float combineDistanceComponents(float componentA, float componentB) const override {
    return componentA + componentB;
  }

  float finalise(float combinedComponents) const override {
    return combinedComponents;
  }
};

}  // namespace ufo

#endif  // UFO_UTILS_GEODESICDISTANCECALCULATOR_H_
