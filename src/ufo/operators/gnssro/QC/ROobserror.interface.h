/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_INTERFACE_H_
#define UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_INTERFACE_H_

#include "ufo/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

typedef int F90roerr;

/// Interface to Fortran RO observation error routines

extern "C" {
  void ufo_roobserror_create_f90(F90roerr &, const ioda::ObsSpace &,
                                 const eckit::Configuration &, const oops::ObsVariables &);
  void ufo_roobserror_delete_f90(F90roerr &);
  void ufo_roobserror_prior_f90(const F90roerr &);
}  // extern C

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_INTERFACE_H_
