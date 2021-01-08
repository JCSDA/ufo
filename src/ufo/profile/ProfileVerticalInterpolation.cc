/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/ProfileVerticalInterpolation.h"

namespace ufo {
  void profileVerticalInterpolation(const std::vector <float> &coordIn,
                                    const std::vector <float> &valuesIn,
                                    const std::vector <float> &coordOut,
                                    std::vector <float> &valuesOut,
                                    const ProfileInterpolation::InterpolationMethod interpMethod,
                                    const ProfileInterpolation::CoordinateOrder coordOrder,
                                    const ProfileInterpolation::OutOfBoundsTreatment outOfBounds)
  {
    const float missingValueFloat = util::missingValue(missingValueFloat);
    // Number of input levels.
    const int nLevsIn = static_cast<int> (coordIn.size());
    // Coordinate multiplier (depends upon the coordinate ordering).
    const float coordMultiplier =
      coordOrder == ProfileInterpolation::CoordinateOrder::Descending ? -1.0 : 1.0;

    // Initial value of 'previous' output coordinate is +/- the maximum float.
    float previousCoordOut = coordMultiplier * std::numeric_limits<float>::max();
    // Initial value of 'previous' input level index is the maximum integer.
    int previousLevIn = std::numeric_limits<int>::max();

    // Loop over each output level.
    for (int jlevOut = 0; jlevOut < coordOut.size(); ++jlevOut) {
      if (coordOut[jlevOut] == missingValueFloat) continue;
      // Index of input level closest to output level (from below or above depending on
      // whether the coordinates are in ascending or descending order).
      int jlevIn = 0;

      // Position at which to start searching for the closest input level.
      // By default this is set to the start of the vector.
      int PosStart = 0;
      // If the output coordinate for this level is larger than the previous value,
      // start searching at the previous level index.
      // (This is in the case of ascending coordinate order; the opposite criterion must be
      // fulfilled for descending order.)
      if (coordMultiplier * coordOut[jlevOut] > coordMultiplier * previousCoordOut)
        PosStart = previousLevIn;

      // Perform the search for the closest input level.
      for (int Pos = PosStart; Pos < nLevsIn; ++Pos) {
        if (coordMultiplier * coordIn[Pos] > coordMultiplier * coordOut[jlevOut]) {
          jlevIn = Pos;
          break;
        }
      }
      // Modify previous input level index and output coordinate.
      previousLevIn = jlevIn;
      previousCoordOut = coordOut[jlevOut];

      // Perform the interpolation.
      if (coordOut[jlevOut] == coordIn[jlevIn]) {
        // Copy the value if the input and output coordinate are identical.
        valuesOut[jlevOut] = valuesIn[jlevIn];
      } else if (outOfBounds != ProfileInterpolation::OutOfBoundsTreatment::Extrapolate &&
                 coordMultiplier * coordOut[jlevOut] > coordMultiplier * coordIn[nLevsIn - 1]) {
        // If the output coordinate is either:
        // - larger than the largest input coordinate if the coordinates are in ascending order, or
        // - smaller than the smallest input coordinate if the coordinates are in descending order
        // and there is no extrapolation, set the output value to one of two possibilities:
        if (outOfBounds == ProfileInterpolation::OutOfBoundsTreatment::SetMissing)
          valuesOut[jlevOut] = missingValueFloat;
        else if (outOfBounds == ProfileInterpolation::OutOfBoundsTreatment::SetToBound)
          valuesOut[jlevOut] = valuesIn[nLevsIn - 1];
      } else if (outOfBounds != ProfileInterpolation::OutOfBoundsTreatment::Extrapolate &&
                 coordMultiplier * coordOut[jlevOut] < coordMultiplier * coordIn[0]) {
        // If the output coordinate is either:
        // - smaller than the smallest input coordinate if the coordinates
        //   are in ascending order, or
        // - larger than the largest input coordinate if the coordinates are in descending order
        // and there is no extrapolation, set the output value to one of two possibilities:
        if (outOfBounds == ProfileInterpolation::OutOfBoundsTreatment::SetMissing)
          valuesOut[jlevOut] = missingValueFloat;
        else if (outOfBounds == ProfileInterpolation::OutOfBoundsTreatment::SetToBound)
          valuesOut[jlevOut] = valuesIn[0];
      } else {
        // Interpolation (or extrapolation if the output coordinate is out-of-bounds)
        // will be performed.
        // Compute the interpolation factor.
        float Interp_factor = 0.0;
        if (interpMethod == ProfileInterpolation::InterpolationMethod::LogLinear) {
          // Log-linear interpolation.
          if (coordOut[jlevOut] > 0.0 &&
              coordIn[jlevIn] > 0.0 &&
              coordIn[jlevIn + 1] > 0.0) {
            Interp_factor = log(coordOut[jlevOut] / coordIn[jlevIn]) /
              log(coordIn[jlevIn + 1] / coordIn[jlevIn]);
          }
        } else {
          // Linear interpolation.
          if (coordIn[jlevIn + 1] - coordIn[jlevIn] != 0) {
            Interp_factor = (coordOut[jlevOut] - coordIn[jlevIn]) /
              (coordIn[jlevIn + 1] - coordIn[jlevIn]);
          }
        }
        // Determine the interpolated value.
        if (valuesIn[jlevIn] != missingValueFloat &&
            valuesIn[jlevIn + 1] != missingValueFloat) {
          valuesOut[jlevOut] = valuesIn[jlevIn] +
            (valuesIn[jlevIn + 1] - valuesIn[jlevIn]) * Interp_factor;
        }
      }
    }
  }
}  // namespace ufo
