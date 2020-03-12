/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ONEDVARFORTRAN_INTERFACE_H_
#define UFO_ONEDVARFORTRAN_INTERFACE_H_

#include "../Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90onedvarcheck;

/// Interface to Fortran routines

extern "C" {
  void ufo_onedvarfortran_create_f90(F90onedvarcheck &, const ioda::ObsSpace &, const eckit::Configuration *, 
                                     const int &, const int &);
  void ufo_onedvarfortran_delete_f90(F90onedvarcheck &);
  void ufo_onedvarfortran_prior_f90(const F90onedvarcheck &, const F90goms &);
  void ufo_onedvarfortran_post_f90(const F90onedvarcheck &, const oops::Variables &, const F90goms &);
}  // extern C

}  // namespace ufo

#endif  // UFO_ONEDVARFORTRAN_INTERFACE_H_
