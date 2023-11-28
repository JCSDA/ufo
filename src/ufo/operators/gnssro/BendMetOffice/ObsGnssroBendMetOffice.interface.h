/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_INTERFACE_H_
#define UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
// Gnssro bending angle observation operators - (Met Office 1D)
// -----------------------------------------------------------------------------
  void ufo_gnssro_bendmetoffice_setup_f90(F90hop &, const bool &, const bool &,
                                          const float &, const int &, const int &, const bool &);
  void ufo_gnssro_bendmetoffice_delete_f90(F90hop &);
  void ufo_gnssro_bendmetoffice_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                       const int &, const int &, double &, const F90goms &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_INTERFACE_H_
