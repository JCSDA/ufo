/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_INTERFACE_H_
#define UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
// Scatwind neutral wind observation operator - (Met Office)
// -----------------------------------------------------------------------------
  void ufo_scatwind_neutralmetoffice_setup_f90(F90hop &,
                                               const bool &,
                                               const int &,
                                               const oops::ObsVariables &,
                                               oops::Variables &,
                                               const size_t &,
                                               const int &,
                                               const eckit::Configuration &);
  void ufo_scatwind_neutralmetoffice_delete_f90(F90hop &);
  void ufo_scatwind_neutralmetoffice_simobs_f90(const F90hop &,
                                                const F90goms &,
                                                const ioda::ObsSpace &,
                                                const int &,
                                                const int &,
                                                double &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_INTERFACE_H_
