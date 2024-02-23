/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "ufo/Fortran.h"

namespace ufo {
namespace fov {

/// Interface to Fortran UFO field of view routines

extern "C" {

// -----------------------------------------------------------------------------

void ufo_fov_setup_f90(F90fov &, const int &, const char *, const int &,
                       const char *, bool &, int &);
void ufo_fov_delete_f90(F90fov &);
void ufo_fov_ellipse_f90(const F90fov &, const int &, const char *,
                         const int &, const double &, const double &,
                         const double &, const int &, double &, double &);
void ufo_antenna_power_within_fov_f90(const F90fov &, const int &, const char *,
                                      const int &, const double &, const double &,
                                      const double &, const double &, const double &,
                                      double &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace fov
}  // namespace ufo
