/*
 * (C) Copyright 2023 NOAA, UCAR, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdexcept>
#include <utility>

#include "ufo/utils/NearestNeighborInterpolation.h"
#include "ufo/utils/VertInterp.interface.h"

namespace ufo {

NearestNeighborInterpolation::NearestNeighborInterpolation(
    std::vector<double> sortedAbscissas, std::vector<double> ordinates) {
  if (sortedAbscissas.empty())
    throw std::invalid_argument("At least one interpolation point must be provided");

  if (sortedAbscissas.size() != ordinates.size())
    throw std::invalid_argument("The number of abscissas must be the same as that of ordinates");

  abscissas_ = std::move(sortedAbscissas);
  ordinates_ = std::move(ordinates);
}

double NearestNeighborInterpolation::operator()(double abscissa) const {
  return interpolate(abscissas_, ordinates_, abscissa);
}

double NearestNeighborInterpolation::interpolate(const std::vector<double> &sortedAbscissas,
                                      const std::vector<double> &ordinates,
                                      double abscissa) {
  if (sortedAbscissas.size() == 1) {
    // The Fortran functions don't handle this case correctly.
    return ordinates[0];
  }

  int idx = 0;
  nearestneighbor_interp_index_f90(sortedAbscissas.size(), abscissa, sortedAbscissas.data(), idx);

  double f = 0.0;
  nearestneighbor_interp_apply_f90(ordinates.size(), ordinates.data(), f, idx);

  return f;
}

}  // namespace ufo
