/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_INTERFACE_H_
#define UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "ufo/Fortran.h"
#include "ufo/Locations.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
// Gnssro bending angle observation operators - (ROPP2D)
// -----------------------------------------------------------------------------
  void ufo_gnssro_2d_locs_init_f90(const F90hop &, F90locs &, const ioda::ObsSpace &);
  void ufo_gnssro_bndropp2d_setup_f90(F90hop &, const eckit::Configuration &, const int &);
  void ufo_gnssro_bndropp2d_delete_f90(F90hop &);
  void ufo_gnssro_bndropp2d_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                       const int &, double &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_INTERFACE_H_
