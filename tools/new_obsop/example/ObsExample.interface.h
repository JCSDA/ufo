/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_EXAMPLE_OBSEXAMPLE_INTERFACE_H_
#define UFO_EXAMPLE_OBSEXAMPLE_INTERFACE_H_

#include "ioda/ObsSpace.h"

#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO example routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_example_setup_f90(F90hop &, const eckit::Configuration * const *,
                             char *, char *, const int &);
  void ufo_example_delete_f90(F90hop &);
  void ufo_example_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const int &, double &, const F90obias &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_EXAMPLE_OBSEXAMPLE_INTERFACE_H_
