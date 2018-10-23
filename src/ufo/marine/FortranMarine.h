/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_MARINE_FORTRANMARINE_H_
#define UFO_MARINE_FORTRANMARINE_H_

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
//  Ice concentration observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_seaicefrac_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicefrac_delete_f90(F90hop &);
  void ufo_seaicefrac_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                 const F90ovec &, const F90obias &);
  void ufo_seaicefrac_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicefrac_tlad_delete_f90(F90hop &);
  void ufo_seaicefrac_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_seaicefrac_simobs_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_seaicefrac_simobs_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ice thickness observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_seaicethick_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicethick_delete_f90(F90hop &);
  void ufo_seaicethick_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                  const F90ovec &, const F90obias &);
  void ufo_seaicethick_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicethick_tlad_delete_f90(F90hop &);
  void ufo_seaicethick_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_seaicethick_simobs_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_seaicethick_simobs_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Steric Height observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_stericheight_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_stericheight_delete_f90(F90hop &);
  void ufo_stericheight_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                   const F90ovec &, const F90obias &);
  void ufo_stericheight_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_stericheight_tlad_delete_f90(F90hop &);
  void ufo_stericheight_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_stericheight_tlad_gettraj_f90(const F90hop &, const int &, const int &, F90goms &);
  void ufo_stericheight_simobs_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_stericheight_simobs_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ocean Insitu Temperature observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_insitutemperature_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_insitutemperature_delete_f90(F90hop &);
  void ufo_insitutemperature_simobs_f90(const F90hop &, const F90goms &,
                                        const ioda::ObsSpace &, const F90ovec &,
                                        const F90obias &);
  void ufo_insitutemperature_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_insitutemperature_tlad_delete_f90(F90hop &);
  void ufo_insitutemperature_tlad_settraj_f90(const F90hop &, const F90goms &,
                                              const ioda::ObsSpace &);
  void ufo_insitutemperature_simobs_tl_f90(const F90hop &, const F90goms &,
                                           const ioda::ObsSpace &, const F90ovec &);
  void ufo_insitutemperature_simobs_ad_f90(const F90hop &, const F90goms &,
                                           const ioda::ObsSpace &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ocean ADT observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_adt_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_adt_delete_f90(F90hop &);
  void ufo_adt_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                          const F90ovec &, const F90obias &);
  void ufo_adt_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_adt_tlad_delete_f90(F90hop &);
  void ufo_adt_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_adt_simobs_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_adt_simobs_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ocean Sea-Surface Temperature observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_seasurfacetemp_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seasurfacetemp_delete_f90(F90hop &);
  void ufo_seasurfacetemp_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                                     const F90ovec &, const F90obias &);
  void ufo_seasurfacetemp_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seasurfacetemp_tlad_delete_f90(F90hop &);
  void ufo_seasurfacetemp_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_seasurfacetemp_simobs_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_seasurfacetemp_simobs_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_MARINE_FORTRANMARINE_H_
