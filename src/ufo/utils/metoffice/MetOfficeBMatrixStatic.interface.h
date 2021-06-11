/*
 * (C) Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_INTERFACE_H_
#define UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_INTERFACE_H_

#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {

/// Interface to Fortran routines
extern "C" {
  void ufo_metoffice_bmatrixstatic_setup_f90(F90obfilter &, const eckit::Configuration &,
                                             size_t &, size_t &);
  void ufo_metoffice_bmatrixstatic_delete_f90(F90obfilter &);
  void ufo_metoffice_bmatrixstatic_getelements_f90(F90obfilter &, const size_t &,
                                                   const size_t &, float *, float *, float *);
}  // extern C

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_INTERFACE_H_
