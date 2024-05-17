/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_INTERFACE_H_
#define TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO example routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_example_tlad_setup_f90(F90hop &, const eckit::Configuration &,
                                  const oops::ObsVariables &, oops::Variables &);
  void ufo_example_tlad_delete_f90(F90hop &);
  void ufo_example_tlad_settraj_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                    const F90goms &);
  void ufo_example_simobs_tl_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                  const int &, double &);
  void ufo_example_simobs_ad_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                  const int &, const double &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_INTERFACE_H_
