/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TOOLS_NEW_QC_EXAMPLE_EXAMPLE_INTERFACE_H_
#define TOOLS_NEW_QC_EXAMPLE_EXAMPLE_INTERFACE_H_

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
  void ufo_example_create_f90(F90check &, const eckit::Configuration &,
                              oops::Variables &);
  void ufo_example_delete_f90(F90check &);
  void ufo_example_prior_f90(const F90check &, const ioda::ObsSpace &,
                             const F90goms &);
  void ufo_example_post_f90(const F90check &, const ioda::ObsSpace &, const int &,
                            const int &, const double &, const double &, const F90goms &);
}  // extern C

}  // namespace ufo

#endif  // TOOLS_NEW_QC_EXAMPLE_EXAMPLE_INTERFACE_H_
