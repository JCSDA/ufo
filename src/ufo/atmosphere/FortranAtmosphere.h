/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ATMOSPHERE_FORTRANATMOSPHERE_H_
#define UFO_ATMOSPHERE_FORTRANATMOSPHERE_H_

#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  Aircraft t observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_aircraft_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aircraft_delete_f90(F90hop &);
  void ufo_aircraft_simobs_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &,
                              const F90obias &);
  void ufo_aircraft_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aircraft_tlad_delete_f90(F90hop &);
  void ufo_aircraft_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_aircraft_simobs_tl_f90(const F90hop &, const F90goms &, const F90odb &,
                                      const F90ovec &);
  void ufo_aircraft_simobs_ad_f90(const F90hop &, const F90goms &, const F90odb &,
                                      const F90ovec &);

// -----------------------------------------------------------------------------
//  Radiosonde t observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_radiosonde_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiosonde_delete_f90(F90hop &);
  void ufo_radiosonde_simobs_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &,
                                const F90obias &);
  void ufo_radiosonde_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiosonde_tlad_delete_f90(F90hop &);
  void ufo_radiosonde_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_radiosonde_simobs_tl_f90(const F90hop &, const F90goms &, const F90odb &,
                                        const F90ovec &);
  void ufo_radiosonde_simobs_ad_f90(const F90hop &, const F90goms &, const F90odb &,
                                        const F90ovec &);

// -----------------------------------------------------------------------------
//  Radiance observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_radiance_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiance_delete_f90(F90hop &);
  void ufo_radiance_simobs_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &,
                            const F90obias &);
  void ufo_radiance_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiance_tlad_delete_f90(F90hop &);
  void ufo_radiance_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_radiance_simobs_tl_f90(const F90hop &, const F90goms &, const F90odb &,
                                    const F90ovec &);
  void ufo_radiance_simobs_ad_f90(const F90hop &, const F90goms &, const F90odb &,
                                    const F90ovec &);

// -----------------------------------------------------------------------------
//  Gnssro observation operators-refractivity (GSI) and their tl/ad
// -----------------------------------------------------------------------------
  void ufo_gnssro_ref_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_gnssro_ref_delete_f90(F90hop &);
  void ufo_gnssro_ref_simobs_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &,
                          const F90obias &);
  void ufo_gnssro_ref_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_gnssro_ref_tlad_delete_f90(F90hop &);
  void ufo_gnssro_ref_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_gnssro_ref_simobs_tl_f90(const F90hop &, const F90goms &, const F90odb &,
                                    const F90ovec &);
  void ufo_gnssro_ref_simobs_ad_f90(const F90hop &, const F90goms &, const F90odb &,
                                    const F90ovec &);

// -----------------------------------------------------------------------------
//  Gnssro observation operators-bending angle (GSI) and their tl/ad
// -----------------------------------------------------------------------------
  void ufo_gnssro_bndgsi_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_gnssro_bndgsi_delete_f90(F90hop &);
  void ufo_gnssro_bndgsi_simobs_f90(const F90hop &, const F90goms &, const F90odb &,
                          const F90ovec &,
                          const F90obias &);
  void ufo_gnssro_bndgsi_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_gnssro_bndgsi_tlad_delete_f90(F90hop &);
  void ufo_gnssro_bndgsi_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_gnssro_bndgsi_simobs_tl_f90(const F90hop &, const F90goms &, const F90odb &,
                                       const F90ovec &);
  void ufo_gnssro_bndgsi_simobs_ad_f90(const F90hop &, const F90goms &, const F90odb &,
                                       const F90ovec &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_ATMOSPHERE_FORTRANATMOSPHERE_H_
