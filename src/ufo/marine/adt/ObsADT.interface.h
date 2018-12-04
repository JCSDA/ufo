/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_MARINE_ADT_OBSADT_INTERFACE_H_
#define UFO_MARINE_ADT_OBSADT_INTERFACE_H_

#include "ioda/ObsSpace.h"

#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO adt routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_adt_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_adt_delete_f90(F90hop &);
  void ufo_adt_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const int &, double &, const F90obias &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_MARINE_ADT_OBSADT_INTERFACE_H_
