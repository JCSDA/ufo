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

// Variables key type
typedef int F90vars;
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

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void ufo_var_create_f90(int & keyVars, const eckit::Configuration * const *);
  void ufo_var_clone_f90(const int & keyVars, int & keyVars_other);
  void ufo_var_info_f90(const int & keyVars, int &);
  void ufo_var_delete_f90(int & keyVars);

// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void ufo_loc_delete_f90(F90locs &);
  void ufo_loc_nobs_f90(const F90locs &, int &);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void ufo_geovals_create_f90(F90goms &);
  void ufo_geovals_delete_f90(F90goms &);
  void ufo_geovals_zero_f90(const F90goms &);
  void ufo_geovals_random_f90(const F90goms &);
  void ufo_geovals_dotprod_f90(const F90goms &, const F90goms &, double &);
  void ufo_geovals_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void ufo_geovals_read_file_f90(const F90goms &, const eckit::Configuration * const *);
  void ufo_geovals_write_file_f90(const F90goms &, const eckit::Configuration * const *);

// -----------------------------------------------------------------------------
//  Wind speed observations
// -----------------------------------------------------------------------------
  void ufo_wspeed_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_wspeed_delete_f90(F90hop &);

  void ufo_wspeed_eqv_f90(const F90goms &, const F90ovec &);
  void ufo_wspeed_equiv_tl_f90(const F90goms &, const F90ovec &, const F90goms &, const double &);
  void ufo_wspeed_equiv_ad_f90(const F90goms &, const F90ovec &, const F90goms &, double &);

  void ufo_wspeed_gettraj_f90(const F90hop &, const int &, F90goms &);
  void ufo_wspeed_settraj_f90(const F90goms &, const F90goms &);
  void ufo_wspeed_inputs_f90(const F90hop &, F90vars &);

// -----------------------------------------------------------------------------
//  Radiance observations
// -----------------------------------------------------------------------------
  void ufo_radiance_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_radiance_delete_f90(F90hop &);

  void ufo_radiance_eqv_f90(const F90goms &, const F90ovec &);
  void ufo_radiance_equiv_tl_f90(const F90goms &, const F90ovec &, const F90goms &, const double &);
  void ufo_radiance_equiv_ad_f90(const F90goms &, const F90ovec &, const F90goms &, double &);

  void ufo_radiance_gettraj_f90(const F90hop &, const int &, F90goms &);  // copied over from wspeed, don't think it's required
  void ufo_radiance_settraj_f90(const F90goms &, const F90goms &);   // copied over from wspeed, don't think it's required
  void ufo_radiance_inputs_f90(const F90hop &, F90vars &);

// -----------------------------------------------------------------------------
//  Conventional q observations
// -----------------------------------------------------------------------------
  void ufo_convq_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_convq_delete_f90(F90hop &);

  void ufo_convq_eqv_f90(const F90goms &, const F90ovec &);
  void ufo_convq_inputs_f90(const F90hop &, F90vars &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void ufo_obsvec_setup_f90(F90ovec &, const F90odb &);
  void ufo_obsvec_clone_f90(const F90ovec &, F90ovec &);
  void ufo_obsvec_delete_f90(F90ovec &);

  void ufo_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void ufo_obsvec_zero_f90(const F90ovec &);
  void ufo_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void ufo_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void ufo_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void ufo_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void ufo_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void ufo_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void ufo_obsvec_invert_f90(const F90ovec &);
  void ufo_obsvec_random_f90(const F90ovec &);
  void ufo_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void ufo_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void ufo_obsvec_nobs_f90(const F90ovec &, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void ufo_obsdb_setup_f90(F90odb &, const eckit::Configuration * const *);
  void ufo_obsdb_delete_f90(F90odb &);
  void ufo_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void ufo_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void ufo_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              F90locs &);
  void ufo_obsdb_getgeovals_f90(const F90odb &, const F90vars &,
                                const util::DateTime * const *, const util::DateTime * const *,
                                F90goms &);
  void ufo_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void ufo_obsdb_nobs_f90(const F90odb &, int &);

// -----------------------------------------------------------------------------
}  // extern C

}  // namespace ufo
#endif  // UFO_FORTRAN_H_
