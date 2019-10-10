/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SURFACE_CORRECTION_HCORRECTION_INTERFACE_H_
#define UFO_SURFACE_CORRECTION_HCORRECTION_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90check;

/// Interface to Fortran routines

extern "C" {
  void ufo_hcorrection_create_f90(F90check &, const eckit::Configuration *,
                              oops::Variables &);
  void ufo_hcorrection_delete_f90(F90check &);
  void ufo_hcorrection_prior_f90(const F90check &, const ioda::ObsSpace &,
                             const F90goms &);
}  // extern C

}  // namespace ufo

#endif  // UFO_SURFACE_CORRECTION_HCORRECTION_INTERFACE_H_
