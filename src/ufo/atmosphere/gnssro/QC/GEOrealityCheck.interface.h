/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_INTERFACE_H_
#define UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_INTERFACE_H_

#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90georealitycheck;

/// Interface to Fortran RO observation geophysical reality check routines

extern "C" {
  void ufo_georealitycheck_create_f90(F90georealitycheck &, const ioda::ObsSpace &, const eckit::Configuration *);
  void ufo_georealitycheck_delete_f90(F90georealitycheck &);
  void ufo_georealitycheck_prior_f90(const F90georealitycheck &, const F90goms &);
  void ufo_georealitycheck_post_f90(const F90georealitycheck &, const int &, const double &);
}  // extern C

}  // namespace ufo

#endif  // UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_INTERFACE_H_
