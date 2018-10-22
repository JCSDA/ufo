/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FORTRANGEOVALS_H_
#define UFO_FORTRANGEOVALS_H_

#include "Fortran.h"

namespace ufo {

/// Interface to Fortran UFO GeoVals routines
/*!
 * The core of the UFO GeoVals is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

  void ufo_geovals_setup_f90(F90goms &, const F90locs &, const eckit::Configuration * const *);
  void ufo_geovals_create_f90(F90goms &);
  void ufo_geovals_delete_f90(F90goms &);
  void ufo_geovals_copy_f90(const F90goms &, F90goms &);
  void ufo_geovals_zero_f90(const F90goms &);
  void ufo_geovals_abs_f90(const F90goms &);
  void ufo_geovals_rms_f90(const F90goms &, double &);
  void ufo_geovals_setup_random_f90(const F90goms &, const eckit::Configuration * const *,
                                    const eckit::Configuration * const *);
  void ufo_geovals_setup_allocone_f90(const F90goms &, const eckit::Configuration * const *,
                                    const eckit::Configuration * const *);
  void ufo_geovals_analytic_init_f90(F90goms &, const F90locs &,
                                     const eckit::Configuration * const *);
  void ufo_geovals_random_f90(const F90goms &);
  void ufo_geovals_scalmult_f90(const F90goms &, const double &);
  void ufo_geovals_assign_f90(const F90goms &, const F90goms &);
  void ufo_geovals_add_f90(const F90goms &, const F90goms &);
  void ufo_geovals_diff_f90(const F90goms &, const F90goms &);
  void ufo_geovals_normalize_f90(const F90goms &, const F90goms &);
  void ufo_geovals_dotprod_f90(const F90goms &, const F90goms &, double &);
  void ufo_geovals_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void ufo_geovals_maxloc_f90(const F90goms &, double &, int &, int &);
  void ufo_geovals_read_file_f90(const F90goms &, const eckit::Configuration * const *,
                                 const eckit::Configuration * const *);
  void ufo_geovals_write_file_f90(const F90goms &, const eckit::Configuration * const *);


}  // extern C

}  // namespace ufo
#endif  // UFO_FORTRANGEOVALS_H_
