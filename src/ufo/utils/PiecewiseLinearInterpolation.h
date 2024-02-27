/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_PIECEWISELINEARINTERPOLATION_H_
#define UFO_UTILS_PIECEWISELINEARINTERPOLATION_H_

#include <tuple>
#include <utility>
#include <vector>

namespace ufo {

/// \brief Represents a piecewise linear interpolation of a set of data points.
class PiecewiseLinearInterpolation
{
 public:
  /// \brief Create an object representing a piecewise linear interpolation of the data points
  /// (sortedAbscissas[i], ordinates[i]).
  ///
  /// Both arguments must have the same length and be non-empty. The elements of \p sortedAbscissas
  /// must be sorted.
  PiecewiseLinearInterpolation(std::vector<double> sortedAbscissas,
                               std::vector<double> ordinates);

  /// \brief Evaluate the interpolated function at \p abscissa.
  double operator()(double abscissa) const;

  /// \brief Convenience function interpolating the data points (sortedAbscissas[i], ordinates[i])
  /// at \p abscissa without creating a PiecewiseLinearInterpolation object.
  static double interpolate(const std::vector<double> &sortedAbscissas,
                            const std::vector<double> &ordinates,
                            double abscissa);

  static std::pair<int, double> interpolationIndexAndWeight
    (const std::vector<double> &sortedAbscissas, double abscissa);

 private:
  std::vector<double> abscissas_;
  std::vector<double> ordinates_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PIECEWISELINEARINTERPOLATION_H_
