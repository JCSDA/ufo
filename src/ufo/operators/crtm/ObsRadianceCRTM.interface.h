/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_OBSRADIANCECRTM_INTERFACE_H_
#define UFO_OPERATORS_CRTM_OBSRADIANCECRTM_INTERFACE_H_

#include "ufo/Fortran.h"

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

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  Radiance observation operator
// -----------------------------------------------------------------------------
  void ufo_radiancecrtm_setup_f90(F90hop &, const eckit::Configuration &,
                                  const int &, const int &,
                                  oops::Variables &, const eckit::mpi::Comm &);
  void ufo_radiancecrtm_delete_f90(F90hop &);
  void ufo_radiancecrtm_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const size_t &, const size_t &, double &, const F90goms &,
                               const void*);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSRADIANCECRTM_INTERFACE_H_
