/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEVERTICALINTERPOLATION_H_
#define UFO_PROFILE_PROFILEVERTICALINTERPOLATION_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"

#include "ufo/utils/Constants.h"

namespace ufo {
  namespace ProfileInterpolation {
    enum class InterpolationMethod {Linear, LogLinear};
    enum class CoordinateOrder {Ascending, Descending};
    enum class OutOfBoundsTreatment {SetMissing, SetToBound, Extrapolate};
  }

  /// Performs vertical interpolation from model levels onto observation levels.
  /// \param[in] coordIn: Coordinates of input values.
  /// \param[in] valuesIn: Input values.
  /// \param[in] coordOut: Coordinates of output values.
  /// \param[out] valuesOut: Interpolated values.
  /// \param[in] interpMethod: Interpolation method (linear or log-linear).
  /// \param[in] coordOrder: Order of coordinates (ascending or desceding).
  /// \param[in] outOfBounds: How to treat out-of-bounds data.
  void profileVerticalInterpolation(const std::vector <float> &coordIn,
                                    const std::vector <float> &valuesIn,
                                    const std::vector <float> &coordOut,
                                    std::vector <float> &valuesOut,
                                    const ProfileInterpolation::InterpolationMethod interpMethod =
                                    ProfileInterpolation::InterpolationMethod::Linear,
                                    const ProfileInterpolation::CoordinateOrder coordOrder =
                                    ProfileInterpolation::CoordinateOrder::Ascending,
                                    const ProfileInterpolation::OutOfBoundsTreatment outOfBounds =
                                    ProfileInterpolation::OutOfBoundsTreatment::SetToBound);
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEVERTICALINTERPOLATION_H_
