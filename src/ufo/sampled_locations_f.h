/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SAMPLED_LOCATIONS_F_H_
#define UFO_SAMPLED_LOCATIONS_F_H_

#include <vector>

#include "ufo/SampledLocations.h"

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to SampledLocations
// -----------------------------------------------------------------------------
namespace ufo {

extern "C" {
  std::size_t sampled_locations_get_npaths_f(const SampledLocations &);
  void sampled_locations_get_lons_f(
      const SampledLocations &, const std::size_t &, double *);
  void sampled_locations_get_lats_f(
      const SampledLocations &, const std::size_t &, double *);
  void sampled_locations_get_timemask_f(
      const SampledLocations &, const util::DateTime &,
      const util::DateTime &, const std::size_t &, bool *);
}

}  // namespace ufo

#endif  // UFO_SAMPLED_LOCATIONS_F_H_
