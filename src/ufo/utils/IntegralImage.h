/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_INTEGRALIMAGE_H_
#define UFO_UTILS_INTEGRALIMAGE_H_

#include "eckit/exception/Exceptions.h"

namespace ufo {

/// Create integral image using a 2D input Eigen array.
template <typename T>
void createIntegralImage(T & array) {
  for (std::size_t i = 1; i < array.rows(); ++i) {
    array(i, 0) += array(i - 1, 0);
  }
  for (std::size_t j = 1; j < array.cols(); ++j) {
    array(0, j) += array(0, j - 1);
  }
  for (std::size_t i = 1; i < array.rows(); ++i) {
    for (std::size_t j = 1; j < array.cols(); ++j) {
      array(i, j) += array(i, j - 1) +
        array(i - 1, j) -
        array(i - 1, j - 1);
    }
  }
}

/// Sum a specified patch of an integral image.
/// No check is made on the bounds of the patch; to do so,
/// use the function `sumIntegralImagePatchWithBoundsCheck`.
template <typename T>
typename T::Scalar sumIntegralImagePatch
  (const T & array,
   const std::size_t i, const std::size_t j,
   const std::size_t istep, const std::size_t jstep) {
  return array(i + istep, j + jstep) -
    array(i + istep, j - jstep - 1) -
    array(i - istep - 1, j + jstep) +
    array(i - istep - 1, j - jstep - 1);
}

/// Sum a specified patch of an integral image.
/// Ensure none of the patch bounds are outside the range
/// of the integral image.
template <typename T>
typename T::Scalar sumIntegralImagePatchWithBoundsCheck
  (const T & array,
   const std::size_t i, const std::size_t j,
   const std::size_t istep, const std::size_t jstep) {
  if (i + istep >= array.rows() ||
      j + jstep >= array.cols() ||
      i - istep - 1 < 0 ||
      j - jstep - 1 < 0) {
    throw eckit::UserError("One or more bounds exceeded", Here());
  }
  return sumIntegralImagePatch(array, i, j, istep, jstep);
}

}  // namespace ufo

#endif  // UFO_UTILS_INTEGRALIMAGE_H_
