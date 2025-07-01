/*
 * (C) Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_INTERFACE_H_
#define UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_INTERFACE_H_

#include "../../Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90onedvarcheck;

/// Interface to Fortran routines

extern "C" {
  void ufo_refractivityonedvarcheck_create_f90(
    F90onedvarcheck &,
    const ioda::ObsSpace &,
    const int &,
    const char *,
    const int &,
    const char *,
    const bool &,
    const float &,
    const float &,
    const float &,
    const float &,
    const int &,
    const float &,
    const bool &,
    const bool &,
    const float &,
    const float &,
    const float &,
    const int &);
  void ufo_refractivityonedvarcheck_delete_f90(F90onedvarcheck &);
  void ufo_refractivityonedvarcheck_apply_f90(
    const F90onedvarcheck &,
    const F90goms &,
    const int &,
    const char &);
}  // extern C

}  // namespace ufo

#endif  // UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_INTERFACE_H_
