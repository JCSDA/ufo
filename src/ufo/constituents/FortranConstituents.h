/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CONSTITUENTS_FORTRANCONSTITUENTS_H_
#define UFO_CONSTITUENTS_FORTRANCONSTITUENTS_H_

#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  AOD observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_aod_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aod_delete_f90(F90hop &);
  void ufo_aod_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &,
                       const F90obias &);
  void ufo_aod_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aod_tlad_delete_f90(F90hop &);
  void ufo_aod_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_aod_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);
  void ufo_aod_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_CONSTITUENTS_FORTRANCONSTITUENTS_H_
