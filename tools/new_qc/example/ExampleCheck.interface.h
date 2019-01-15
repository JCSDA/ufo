/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_EXAMPLECHECK_INTERFACE_H_
#define UFO_EXAMPLECHECK_INTERFACE_H_

#include "Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90check;

/// Interface to Fortran routines

extern "C" {
  void ufo_examplecheck_create_f90(F90check &, const ioda::ObsSpace &, const eckit::Configuration *);
  void ufo_examplecheck_delete_f90(F90check &);
  void ufo_examplecheck_prior_f90(const F90check &, const F90goms &);
  void ufo_examplecheck_post_f90(const F90check &, const int &, const double &);
}  // extern C

}  // namespace ufo

#endif  // UFO_EXAMPLECHECK_INTERFACE_H_
