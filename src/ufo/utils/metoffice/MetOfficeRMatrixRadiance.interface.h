/*
 * (C) Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_INTERFACE_H_
#define UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_INTERFACE_H_

#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {

// Interface to Fortran routines
extern "C" {
  void ufo_metoffice_rmatrixradiance_setup_f90(const F90obfilter &, const eckit::Configuration &,
                                               size_t &, size_t &, size_t &);
  void ufo_metoffice_rmatrixradiance_delete_f90(const F90obfilter &);
  void ufo_metoffice_rmatrixradiance_getelements_f90(const F90obfilter &,
                                             const size_t &, int *, float *);
  void ufo_metoffice_rmatrixradiance_print_f90(const F90obfilter &);
}  // extern C

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_INTERFACE_H_
