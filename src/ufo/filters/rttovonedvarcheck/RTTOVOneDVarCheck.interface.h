/*
 * (C) Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_INTERFACE_H_
#define UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_INTERFACE_H_

#include "ufo/Fortran.h"

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
  void ufo_rttovonedvarcheck_create_f90(F90onedvarcheck &, const ioda::ObsSpace &,
                      const eckit::Configuration *, const int &, const int &,
                      oops::Variables &, const int &, const int &);
  void ufo_rttovonedvarcheck_delete_f90(F90onedvarcheck &);
  void ufo_rttovonedvarcheck_apply_f90(const F90onedvarcheck &, const oops::Variables &,
                                       const oops::Variables &, const F90goms &, const int &,
                                       const char &);
}  // extern C

}  // namespace ufo

#endif  // UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_INTERFACE_H_
