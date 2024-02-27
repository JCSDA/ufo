/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdexcept>
#include <utility>

#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/utils/VertInterp.interface.h"

namespace ufo {

PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double> sortedAbscissas, std::vector<double> ordinates) {
  if (sortedAbscissas.empty())
    throw std::invalid_argument("At least one interpolation point must be provided");

  if (sortedAbscissas.size() != ordinates.size())
    throw std::invalid_argument("The number of abscissas must be the same as that of ordinates");

  abscissas_ = std::move(sortedAbscissas);
  ordinates_ = std::move(ordinates);
}

double PiecewiseLinearInterpolation::operator()(double abscissa) const {
  return interpolate(abscissas_, ordinates_, abscissa);
}

double PiecewiseLinearInterpolation::interpolate(const std::vector<double> &sortedAbscissas,
                                                 const std::vector<double> &ordinates,
                                                 double abscissa) {
  if (sortedAbscissas.size() == 1) {
    // The Fortran functions don't handle this case correctly.
    return ordinates[0];
  }

  int wi = 0;
  double wf = 0.0;
  vert_interp_weights_f90(sortedAbscissas.size(), abscissa, sortedAbscissas.data(), wi, wf);

  double f = 0.0;
  vert_interp_apply_f90(ordinates.size(), ordinates.data(), f, wi, wf);

  return f;
}

std::pair<int, double> PiecewiseLinearInterpolation::interpolationIndexAndWeight
(const std::vector<double> &sortedAbscissas, double abscissa) {
  if (sortedAbscissas.size() == 1) {
    // The Fortran functions don't handle this case correctly.
    return {0, 1.0};
  }

  int wi = 0;
  double wf = 0.0;
  vert_interp_weights_f90(sortedAbscissas.size(), abscissa, sortedAbscissas.data(), wi, wf);

  // Convert Fortran index to C++ index
  return {wi - 1, wf};
}
}  // namespace ufo
