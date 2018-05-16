/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FORTRAN_H_
#define UFO_FORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace ufo {

// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Observation vector key type
typedef int F90ovec;
// Obs operator key type
typedef int F90hop;
// Observation space type
typedef int F90odb;
// Observation check key type
typedef int F90ocheck;
// Observation bias key type
typedef int F90obias;

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void ufo_geovals_setup_f90(F90goms &, const F90locs &, const eckit::Configuration * const *);
  void ufo_geovals_create_f90(F90goms &);
  void ufo_geovals_delete_f90(F90goms &);
  void ufo_geovals_zero_f90(const F90goms &);
  void ufo_geovals_setup_random_f90(const F90goms &, const eckit::Configuration * const *,
                                    const eckit::Configuration * const *);
  void ufo_geovals_random_f90(const F90goms &);
  void ufo_geovals_scalmult_f90(const F90goms &, const double &);
  void ufo_geovals_assign_f90(const F90goms &, const F90goms &);
  void ufo_geovals_add_f90(const F90goms &, const F90goms &);    
  void ufo_geovals_dotprod_f90(const F90goms &, const F90goms &, double &);
  void ufo_geovals_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void ufo_geovals_read_file_f90(const F90goms &, const eckit::Configuration * const *,
                                 const eckit::Configuration * const *);
  void ufo_geovals_write_file_f90(const F90goms &, const eckit::Configuration * const *);

// -----------------------------------------------------------------------------
//  Check Observations
// -----------------------------------------------------------------------------
  void ufo_obscheck_setup_f90(F90ocheck &, const eckit::Configuration * const *);
  void ufo_obscheck_delete_f90(F90ocheck &);
  void ufo_postFilter_f90(const F90goms &, const F90ovec &, const F90odb &);
  void ufo_priorFilter_f90(const F90odb &);

// -----------------------------------------------------------------------------
//  Radiosonde t observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_radiosonde_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiosonde_delete_f90(F90hop &);
  void ufo_radiosonde_t_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_radiosonde_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiosonde_tlad_delete_f90(F90hop &);
  void ufo_radiosonde_tlad_settraj_f90(const F90hop &, const F90goms &, const F90odb &);
  void ufo_radiosonde_tlad_t_eqv_tl_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);
  void ufo_radiosonde_tlad_t_eqv_ad_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Radiance observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_radiance_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiance_delete_f90(F90hop &);
  void ufo_radiance_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_radiance_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiance_tlad_delete_f90(F90hop &);
  void ufo_radiance_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_radiance_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);
  void ufo_radiance_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ice concentration observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_seaicefrac_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicefrac_delete_f90(F90hop &);
  void ufo_seaicefrac_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_seaicefrac_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicefrac_tlad_delete_f90(F90hop &);
  void ufo_seaicefrac_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_seaicefrac_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_seaicefrac_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Ice thickness observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_seaicethick_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicethick_delete_f90(F90hop &);
  void ufo_seaicethick_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_seaicethick_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_seaicethick_tlad_delete_f90(F90hop &);
  void ufo_seaicethick_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_seaicethick_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_seaicethick_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90ovec &);

// -----------------------------------------------------------------------------
//  Steric Height observation operator and its tl/ad
// -----------------------------------------------------------------------------
  void ufo_stericheight_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_stericheight_delete_f90(F90hop &);
  void ufo_stericheight_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_stericheight_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_stericheight_tlad_delete_f90(F90hop &);
  void ufo_stericheight_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_stericheight_tlad_gettraj_f90(const F90hop &, const int &, const int &, F90goms &);
  void ufo_stericheight_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90ovec &);
  void ufo_stericheight_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90ovec &);
 
// -----------------------------------------------------------------------------
//  AOD observation operator and its tl/ad 
// -----------------------------------------------------------------------------
  void ufo_aod_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aod_delete_f90(F90hop &);
  void ufo_aod_eqv_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &, const F90obias &);
  void ufo_aod_tlad_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_aod_tlad_delete_f90(F90hop &);
  void ufo_aod_tlad_settraj_f90(const F90hop &, const F90goms &);
  void ufo_aod_tlad_eqv_tl_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);
  void ufo_aod_tlad_eqv_ad_f90(const F90hop &, const F90goms &, const F90odb &, const F90ovec &);

// -----------------------------------------------------------------------------  

}  // extern C

}  // namespace ufo
#endif  // UFO_FORTRAN_H_
