/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_MAXNORMDISTANCECALCULATOR_H_
#define UFO_UTILS_MAXNORMDISTANCECALCULATOR_H_

#include <algorithm>
#include <cmath>

#include "ufo/utils/DistanceCalculator.h"

namespace ufo {

/// A DistanceCalculator implementing the distance function that maps each pair of points
/// (x1, x2, x3, ...) and (y1, y2, y3, ...) to max_i(|x_i - y_i|/w_i), where w_i (i >= 1)
/// are scaling factors.
class MaxNormDistanceCalculator : public DistanceCalculator
{
 public:
  float spatialDistanceComponent(float obsLatitude, float obsLongitude,
                                 float latitudeBinCenter, float longitudeBinCenter,
                                 float inverseLatitudeBinWidth,
                                 float inverseLongitudeBinWidth) const override {
    float latitudeComponent =
        std::abs(obsLatitude - latitudeBinCenter) * inverseLatitudeBinWidth;
    float longitudeComponent =
        std::abs(obsLongitude - longitudeBinCenter) * inverseLongitudeBinWidth;
    return combineDistanceComponents(latitudeComponent, longitudeComponent);
  }

  float nonspatialDistanceComponent(float obs, float binCenter,
                                    float inverseBinWidth) const override {
    return std::abs(obs - binCenter) * inverseBinWidth;
  }

  float combineDistanceComponents(float componentA, float componentB) const override {
    return std::max(componentA, componentB);
  }

  float finalise(float combinedComponents) const override {
    return combinedComponents;
  }
};

}  // namespace ufo

#endif  // UFO_UTILS_MAXNORMDISTANCECALCULATOR_H_
