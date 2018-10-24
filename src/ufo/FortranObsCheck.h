/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FORTRANOBSCHECK_H_
#define UFO_FORTRANOBSCHECK_H_

#include "Fortran.h"

#include "ioda/ObsSpace.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  Check Observations
// -----------------------------------------------------------------------------
  void ufo_obscheck_setup_f90(F90ocheck &, const eckit::Configuration * const *);
  void ufo_obscheck_delete_f90(F90ocheck &);
  void ufo_postFilter_f90(const F90goms &, const int &, const double &, const ioda::ObsSpace &);
  void ufo_priorFilter_f90(const ioda::ObsSpace &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_FORTRANOBSCHECK_H_
