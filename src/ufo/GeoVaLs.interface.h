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

namespace util {
  template <typename T> struct Range;
}

namespace ufo {
  class SampledLocations;

/// Interface to Fortran UFO GeoVals routines
/*!
 * The core of the UFO GeoVals is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
  void ufo_geovals_default_constr_f90(F90goms &);

  /// \brief Creates Fortran GeoVaLs without allocating arrays storing GeoVaL values.
  ///
  /// Deprecated; use ufo_geovals_setup_f90 instead.
  ///
  /// \param key
  ///   Key of the GeoVaLs object to act on.
  /// \param nlocs
  ///   Number of observation locations.
  /// \param vars
  ///   Names of the variables to be stored in the new GeoVaLs in the sampled format.
  /// \param nvars
  ///   Number of variables in `vars`.
  /// \param nsampling_methods
  ///   Number of distinct methods of sampling the observation locations with sets of interpolation
  ///   paths.
  /// \param sampling_method_by_var
  ///   An array of length `nvars` whose ith element is the index of the observation location
  ///   sampling method producing the set of paths along which the ith variable in `vars` will be
  ///   interpolated. Valid values are integers from 0 to `nsampling_methods - 1`.
  /// \param reduced_vars
  ///   Names of the variables to be stored in the new GeoVaLs in the reduced format.
  /// \param is_sampling_method_trivial
  ///   Array of length `nsampling_methods` whose ith element indicates whether the ith sampling
  ///   method produces a set of interpolation paths sampling all locations exactly once in
  ///   ascending order.
  ///
  /// This call must be followed by calls to ufo_geovals_setup_(trivial_)sampling_method_f90 mapping
  /// locations to interpolation paths and ufo_geovals_allocate_f90 allocating specific GeoVaLs.
  void ufo_geovals_partial_setup_f90(F90goms & key, const size_t & nlocs,
                                     const oops::Variables & vars,
                                     const size_t & nvars,
                                     const size_t & nsampling_methods,
                                     const size_t * sampling_method_by_var,
                                     const oops::Variables & reduced_vars,
                                     const bool * is_sampling_method_trivial);

  /// \brief Creates and allocates Fortran GeoVaLs.
  ///
  /// \param key
  ///   Key of the GeoVaLs object to act on.
  /// \param nlocs
  ///   Number of observation locations.
  /// \param vars
  ///   Names of the variables to be stored in the new GeoVaLs in the sampled format.
  /// \param nvars
  ///   Number of variables in `vars`.
  /// \param nlevs
  ///   Array of length `nvars` whose ith element indicates how many values per interpolation path
  ///   will be stored in the GeoVaL corresponding to the ith variable.
  /// \param nsampling_methods
  ///   Number of distinct methods of sampling the observation locations with sets of interpolation
  ///   paths.
  /// \param npaths_by_method
  ///   Array of length `nsampling_methods` whose ith element is the number of paths in the set
  ///   of interpolation paths produced by the ith sampling method.
  /// \param sampling_method_by_var
  ///   Array of length `nvars` whose ith element is the index of the observation location
  ///   sampling method producing the set of paths along which the ith variable will be
  ///   interpolated. Valid values are integers from 0 to `nsampling_methods - 1`.
  /// \param reduced_vars
  ///   Names of the variables to be stored in the new GeoVaLs in the reduced format.
  /// \param nreduced_vars
  ///   Number of variables in `reduced_vars`.
  /// \param nreduced_levs
  ///   Array of length `nreduced_vars` whose ith element indicates how many values per
  ///   location will be stored in the GeoVaL corresponding to the ith variable.
  /// \param is_sampling_method_trivial
  ///   Array of length `nsampling_methods` whose ith element indicates whether the ith sampling
  ///   method produces a set of interpolation paths sampling all locations exactly once in
  ///   ascending order.
  ///
  /// This call must be followed by calls to ufo_geovals_setup_(trivial_)sampling_method_f90 mapping
  /// locations to interpolation paths.
  void ufo_geovals_setup_f90(F90goms & key, const size_t & nlocs,
                             const oops::Variables & vars,
                             const size_t & nvars, const size_t * nlevs,
                             const size_t & nsampling_methods,
                             const size_t * npaths_by_method,
                             const size_t * sampling_method_by_var,
                             const oops::Variables & reduced_vars,
                             const size_t & nreduced_vars, const size_t * nreduced_levs,
                             const bool * is_sampling_method_trivial);

  /// Deprecated, rely on ufo_geovals_setup_f90 to allocate GeoVaLs instead.
  /// Allocates GeoVaLs for \p vars variables with \p nlevels number of levels.
  /// If the GeoVaLs for this variable were allocated before with different size,
  /// aborts.
  void ufo_geovals_allocate_f90(const F90goms &, const size_t & nlevels,
                                const oops::Variables & vars);

  /// \brief Specify which interpolation paths produced by a given method sample which observation
  /// locations.
  ///
  /// \param key
  ///   Key of the GeoVaLs object to act on.
  /// \param sampling_method
  ///   0-based index of the observation location sampling method to set up.
  /// \param npaths
  ///   Number of interpolation paths produced by this sampling method.
  /// \param nlocs
  ///   Number of observation locations.
  /// \param paths_by_loc
  ///   An array of length `nlocs` mapping the index of each observation location to the range of
  ///   (0-based) indices of the paths sampling that location. Specifically, the ith location is
  ///   deemed to be sampled by the paths with indices ranging from `paths_by_loc[i].begin` up to
  ///   but not including `paths_by_loc[i].end`.
  void ufo_geovals_setup_sampling_method_f90(
      F90goms & key, const size_t & sampling_method,
      const size_t & npaths, const size_t & nlocs, const util::Range<size_t> * paths_by_loc);

  /// \brief Designate an observation location sampling method as "trivial", i.e. one producing a
  /// set of interpolation paths such that each location is sampled solely by the path with the
  /// same index.
  ///
  /// \param key
  ///   Key of the GeoVaLs object to act on.
  /// \param sampling_method
  ///   0-based index of the observation location sampling method to designate as "trivial".
  void ufo_geovals_setup_trivial_sampling_method_f90(F90goms & key, const size_t & sampling_method);

  /// \brief Instruct a GeoVaLs object to store some GeoVaLs in the reduced format.
  ///
  /// Space will be allocated for the reduced representation of these variables unless they are
  /// already stored in the sampled format and their reduced and sampled representations are
  /// identical.
  ///
  /// \param key
  ///   Key of the GeoVaLs object to act on.
  /// \param vars
  ///   Names of the variables to be stored in the reduced format. They must not already be stored
  ///   in this format.
  /// \param nvars
  ///   Number of variables in `vars`.
  /// \param nlevs
  ///   Array of length `nvars` whose ith element indicates how many values per
  ///   location will be stored in the GeoVaL corresponding to the ith variable in `vars`.
  void ufo_geovals_add_reduced_vars_f90(const F90goms & key, const oops::Variables & vars,
                                        const size_t & nvars, const size_t * nlevs);

  void ufo_geovals_get_vars_f90(const F90goms &, oops::Variables &, const int &);

  int ufo_geovals_get_default_format_f90(const F90goms &);
  void ufo_geovals_set_default_format_f90(F90goms &, const int &);
  bool ufo_geovals_are_reduced_and_sampled_formats_aliased_f90(const F90goms &,
                                                                const int &, const char *);

  void ufo_geovals_delete_f90(F90goms &);
  void ufo_geovals_copy_f90(const F90goms &, F90goms &);
  void ufo_geovals_copy_one_f90(F90goms &, const F90goms &, const int &);
  void ufo_geovals_zero_f90(const F90goms &);
  void ufo_geovals_reorderzdir_f90(const F90goms &, const int &, const char *,
                                   const int &, const char *);
  void ufo_geovals_abs_f90(const F90goms &);
  void ufo_geovals_rms_f90(const F90goms &, double &);
  void ufo_geovals_analytic_init_f90(F90goms &, const SampledLocations &,
                                     const eckit::Configuration &);
  void ufo_geovals_random_f90(const F90goms &);
  void ufo_geovals_scalmult_f90(const F90goms &, const double &);
  void ufo_geovals_profmult_f90(const F90goms &, const int &, const float &);
  void ufo_geovals_assign_f90(const F90goms &, const F90goms &);
  void ufo_geovals_add_f90(const F90goms &, const F90goms &);
  void ufo_geovals_diff_f90(const F90goms &, const F90goms &);
  void ufo_geovals_schurmult_f90(const F90goms &, const F90goms &);
  void ufo_geovals_normalize_f90(const F90goms &, const F90goms &);
  void ufo_geovals_minmaxavg_f90(const F90goms &, int &, int &, double &, double &, double &);
  void ufo_geovals_maxloc_f90(const F90goms &, double &, int &, int &);
  void ufo_geovals_nlocs_f90(const F90goms &, size_t &);
  void ufo_geovals_nprofiles_f90(const F90goms &, const int &, const char *, const int &, size_t &);
  void ufo_geovals_nlevs_f90(const F90goms &, const int &, const char *, const int &, int &);
  void ufo_geovals_get2d_f90(const F90goms &, const int &, const char *, const int &, const int &,
                             double &);
  void ufo_geovals_get_f90(const F90goms &, const int &, const char *, const int &, const int &,
                           const int &, float &);
  void ufo_geovals_get_profile_f90(const F90goms &, const int &, const char *, const int &,
                                   const int &, const int &, double &);
  void ufo_geovals_getdouble_f90(const F90goms &, const int &, const char *, const int &,
                                 const int &, const int &, double &);
  void ufo_geovals_putdouble_f90(const F90goms &, const int &, const char *, const int &,
                                 const int &, const int &, const double &);
  void ufo_geovals_put_profile_f90(const F90goms &, const int &, const char *, const int &,
                                   const int &, const int &, const double &);
  void ufo_geovals_get_profile_indices_grouped_by_loc_f90(
      const F90goms &key, const int &lvar, const char *var, const int &format,
      const size_t & nlocs, const util::Range<size_t> *profile_indices_grouped_by_location);
  void ufo_geovals_read_file_f90(const F90goms &,
                                 const eckit::Configuration &, const ioda::ObsSpace &,
                                 const oops::Variables &);
  void ufo_geovals_write_file_f90(const F90goms &, const eckit::Configuration &, const size_t &);
  void ufo_geovals_fill_f90(const int &, const int &, const char *, const int &, const int *,
                            const int &, const double *, const bool &);
  void ufo_geovals_fillad_f90(const int &, const int &, const char *, const int &, const int *,
                              const int &, double *, const bool &);
}  // extern C

}  // namespace ufo
#endif  // UFO_GEOVALS_INTERFACE_H_
