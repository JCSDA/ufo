/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GNSSRO_QC_ROGEOREALITYCHECK_INTERFACE_H_
#define UFO_GNSSRO_QC_ROGEOREALITYCHECK_INTERFACE_H_

#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90rogeorealitycheck;

/// Interface to Fortran RO observation geophysical reality check routines

extern "C" {
  void ufo_rogeorealitycheck_create_f90(F90rogeorealitycheck &, const ioda::ObsSpace &,
                                        const eckit::Configuration *);
  void ufo_rogeorealitycheck_delete_f90(F90rogeorealitycheck &);
  void ufo_rogeorealitycheck_prior_f90(const F90rogeorealitycheck &, const F90goms &);
  void ufo_rogeorealitycheck_post_f90(const F90rogeorealitycheck &, const int &, const double &);
}  // extern C

}  // namespace ufo

#endif  // UFO_GNSSRO_QC_ROGEOREALITYCHECK_INTERFACE_H_
