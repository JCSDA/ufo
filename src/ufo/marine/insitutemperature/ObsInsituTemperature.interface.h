/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_
#define UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_

#include "ioda/ObsSpace.h"

#include "ufo/Fortran.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ufo {

/// Interface to Fortran UFO marine/insitutemperature routines

extern "C" {

// -----------------------------------------------------------------------------

  void ufo_insitutemperature_setup_f90(F90hop &, const eckit::Configuration * const *);
  void ufo_insitutemperature_delete_f90(F90hop &);
  void ufo_insitutemperature_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const int &, double &, const F90obias &);
  void ufo_insitutemperature_locateobs_f90(const F90hop &, const ioda::ObsSpace &,
                                           const util::DateTime * const *,
                                           const util::DateTime * const *, F90locs &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_
