/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GEOVALS_INTERFACE_H_
#define UFO_GEOVALS_INTERFACE_H_

#include "Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace oops {
  class Variables;
}


namespace ufo {
  class Locations;

/// Interface to Fortran UFO GeoVals routines
/*!
 * The core of the UFO GeoVals is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
  void ufo_geovals_default_constr_f90(F90goms &);
  /// Deprecated, use ufo_geovals_setup_f90 instead.
  /// Creates Fortran GeoVaLs with key \p key, \p nlocs number of locations
  /// for \p vars variables. Does not allocate internal structures.
  void ufo_geovals_partial_setup_f90(F90goms & key, const size_t & nlocs,
                                     const oops::Variables & vars);
  /// Creates and allocates Fortran GeoVaLs with key \p key, \p nlocs number of
  /// locations for \p vars variables. \p nlevs is a pointer to the first element
  /// of array size \p nvars that contains number of values to allocate for each of
  /// the variables.
  void ufo_geovals_setup_f90(F90goms & key, const size_t & nlocs,
                             const oops::Variables & vars,
                             const size_t & nvars, const size_t & nlevs);
  /// Deprecated, rely on ufo_geovals_setup_f90 to allocate GeoVaLs instead.
  /// Allocates GeoVaLs for \p vars variables with \p nlevels number of levels.
  /// If the GeoVaLs for this variable were allocated before with different size,
  /// aborts.
  void ufo_geovals_allocate_f90(const F90goms &, const size_t & nlevels,
                                const oops::Variables & vars);
  void ufo_geovals_delete_f90(F90goms &);
  void ufo_geovals_copy_f90(const F90goms &, F90goms &);
  void ufo_geovals_copy_one_f90(F90goms &, const F90goms &, const int &);
  void ufo_geovals_zero_f90(const F90goms &);
  void ufo_geovals_reorderzdir_f90(const F90goms &, const int &, const char *,
                                   const int &, const char *);
  void ufo_geovals_abs_f90(const F90goms &);
  void ufo_geovals_rms_f90(const F90goms &, double &);
  void ufo_geovals_analytic_init_f90(F90goms &, const Locations &,
                                     const eckit::Configuration &);
  void ufo_geovals_random_f90(const F90goms &);
  void ufo_geovals_scalmult_f90(const F90goms &, const double &);
  void ufo_geovals_profmult_f90(const F90goms &, const int &, const float &);
  void ufo_geovals_assign_f90(const F90goms &, const F90goms &);
  void ufo_geovals_add_f90(const F90goms &, const F90goms &);
  void ufo_geovals_diff_f90(const F90goms &, const F90goms &);
  void ufo_geovals_schurmult_f90(const F90goms &, const F90goms &);
  void ufo_geovals_normalize_f90(const F90goms &, const F90goms &);
  void ufo_geovals_split_f90(const F90goms &, const F90goms &, const F90goms &);
  void ufo_geovals_merge_f90(const F90goms &, const F90goms &, const F90goms &);
  void ufo_geovals_minmaxavg_f90(const F90goms &, int &, int &, double &, double &, double &);
  void ufo_geovals_maxloc_f90(const F90goms &, double &, int &, int &);
  void ufo_geovals_nlocs_f90(const F90goms &, size_t &);
  void ufo_geovals_nlevs_f90(const F90goms &, const int &, const char *, int &);
  void ufo_geovals_get2d_f90(const F90goms &, const int &, const char *, const int &,
                           double &);
  void ufo_geovals_get_f90(const F90goms &, const int &, const char *, const int &,
                           const int &, float &);
  void ufo_geovals_get_loc_f90(const F90goms &, const int &, const char *, const int &,
                               const int &, double &);
  void ufo_geovals_getdouble_f90(const F90goms &, const int &, const char *, const int &,
                                 const int &, double &);
  void ufo_geovals_putdouble_f90(const F90goms &, const int &, const char *, const int &,
                                 const int &, const double &);
  void ufo_geovals_put_loc_f90(const F90goms &, const int &, const char *, const int &,
                               const int &, const double &);
  void ufo_geovals_read_file_f90(const F90goms &,
                                 const eckit::Configuration &,
                                 const ioda::ObsSpace &, const oops::Variables &);
  void ufo_geovals_write_file_f90(const F90goms &, const eckit::Configuration &, const size_t &);
  void ufo_geovals_fill_f90(const int &, const int &, const int &,
                            const int &, const double &, const bool &);
  void ufo_geovals_fillad_f90(const int &, const int &, const int &,
                              const int &, double &, const bool &);
}  // extern C

}  // namespace ufo
#endif  // UFO_GEOVALS_INTERFACE_H_
