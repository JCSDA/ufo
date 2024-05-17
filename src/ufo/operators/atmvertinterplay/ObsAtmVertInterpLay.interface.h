/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_INTERFACE_H_
#define UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO atmvertinterplay routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_atmvertinterplay_setup_f90(F90hop &, const eckit::Configuration &,
                             const oops::ObsVariables &, oops::Variables &);
  void ufo_atmvertinterplay_delete_f90(F90hop &);
  void ufo_atmvertinterplay_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                             const int &, const int &, double &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_INTERFACE_H_
