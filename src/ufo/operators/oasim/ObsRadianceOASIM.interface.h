/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_INTERFACE_H_
#define UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_INTERFACE_H_

#include "ioda/ObsSpace.h"

#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran OASIM routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_radianceoasim_setup_f90(F90hop &, const eckit::Configuration &,
                                  const int &, const int &);
  void ufo_radianceoasim_delete_f90(F90hop &);
  void ufo_radianceoasim_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                    const int &, const int&, double &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_INTERFACE_H_

