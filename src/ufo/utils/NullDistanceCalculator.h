/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_NULLDISTANCECALCULATOR_H_
#define UFO_UTILS_NULLDISTANCECALCULATOR_H_

#include "ufo/utils/DistanceCalculator.h"

namespace ufo {

/// A DistanceCalculator subclass for cases when you want to go through all the Gaussian Thinning
///  groupObservationsBy... functions in order to group obs into bins, but without calculating any
///  distance norm - e.g. for the select_median option, which does not need distance norm.
class NullDistanceCalculator : public DistanceCalculator {
 public:
  float spatialDistanceComponent(float obsLatitude, float obsLongitude,
                                 float latitudeBinCenter, float longitudeBinCenter,
                                 float inverseLatitudeBinWidth,
                                 float inverseLongitudeBinWidth) const override {
    return 0.0;
  }

  float nonspatialDistanceComponent(float obs, float binCenter,
                                    float inverseBinWidth) const override {
    return 0.0;
  }

  float combineDistanceComponents(float componentA, float componentB) const override {
    return 0.0;
  }

  float finalise(float combinedComponents) const override {
    return 0.0;
  }
};

}  // namespace ufo

#endif  // UFO_UTILS_NULLDISTANCECALCULATOR_H_
