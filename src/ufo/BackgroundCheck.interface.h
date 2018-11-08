/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_BACKGROUNDCHECK_INTERFACE_H_
#define UFO_BACKGROUNDCHECK_INTERFACE_H_

#include "Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90bgcheck;

/// Interface to Fortran background check routines

extern "C" {
  void ufo_bgcheck_create_f90(F90bgcheck &, const ioda::ObsSpace &, const eckit::Configuration *);
  void ufo_bgcheck_delete_f90(F90bgcheck &);
  void ufo_bgcheck_prior_f90(const F90bgcheck &, const F90goms &);
  void ufo_bgcheck_post_f90(const F90bgcheck &, const int &, const double &);
}  // extern C

}  // namespace ufo

#endif  // UFO_BACKGROUNDCHECK_INTERFACE_H_
