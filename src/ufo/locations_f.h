/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LOCATIONS_F_H_
#define UFO_LOCATIONS_F_H_

#include <vector>

#include "ufo/Locations.h"

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to Locations
// -----------------------------------------------------------------------------
namespace ufo {

extern "C" {
  std::size_t locations_get_nlocs_f(const Locations &);
  void locations_get_lons_f(const Locations &, const std::size_t &, double *);
  void locations_get_lats_f(const Locations &, const std::size_t &, double *);
  void locations_get_timemask_f(const Locations &, const util::DateTime &,
                                const util::DateTime &, const std::size_t &, bool *);
}

}  // namespace ufo

#endif  // UFO_LOCATIONS_F_H_
