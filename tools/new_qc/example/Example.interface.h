/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_EXAMPLE_INTERFACE_H_
#define UFO_FILTERS_EXAMPLE_INTERFACE_H_

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
  void ufo_example_create_f90(F90check &, const eckit::Configuration *);
  void ufo_example_delete_f90(F90check &);
  void ufo_example_prior_f90(const F90check &, const ioda::ObsSpace &, const F90goms &);
  void ufo_example_post_f90(const F90check &, const ioda::ObsSpace &, const int &, const double &);
}  // extern C

}  // namespace ufo

#endif  // UFO_FILTERS_EXAMPLE_INTERFACE_H_
