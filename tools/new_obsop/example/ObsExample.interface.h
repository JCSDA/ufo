/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_INTERFACE_H_
#define TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_INTERFACE_H_

#include "ufo/Fortran.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace oops {
  class ObsVariables;
  class Variables;
}

namespace ufo {

/// Interface to Fortran UFO example routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_example_setup_f90(F90hop &, const eckit::Configuration &,
                             const oops::ObsVariables &, oops::Variables &);
  void ufo_example_delete_f90(F90hop &);
  void ufo_example_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const int &, const int &, double &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_INTERFACE_H_
