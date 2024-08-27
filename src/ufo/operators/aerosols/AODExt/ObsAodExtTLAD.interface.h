/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_INTERFACE_H_
#define UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO aodext routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_aodext_tlad_setup_f90(F90hop &, const eckit::Configuration &,
                                  const int &, const int &,
                                  const oops::ObsVariables &, oops::Variables &);
  void ufo_aodext_tlad_delete_f90(F90hop &);
  void ufo_aodext_tlad_settraj_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &);
  void ufo_aodext_simobs_tl_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                  const int &, const int &, double &);
  void ufo_aodext_simobs_ad_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                  const int &, const int &, const double &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_INTERFACE_H_
