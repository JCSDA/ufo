/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_TIMEOPER_OBSTIMEOPER_INTERFACE_H_
#define UFO_TIMEOPER_OBSTIMEOPER_INTERFACE_H_

#include "ioda/ObsSpace.h"

#include "ufo/Fortran.h"
#include "ufo/Locations.h"

namespace ufo {

/// Interface to Fortran UFO timeoper routines

extern "C" {

// -----------------------------------------------------------------------------
  void ufo_timeoper_locs_init_f90(const F90hop &, F90locs &,
                                  const ioda::ObsSpace &,
                                  const util::DateTime * const *,
                                  const util::DateTime * const *,
                                  const util::DateTime * const *,
                                  const util::DateTime * const *,
                                  const util::DateTime * const *);

  void ufo_timeoper_setup_f90(F90hop &, const eckit::Configuration * const *,
                              const int &, const int &);
  void ufo_timeoper_set_timeweight_f90(const F90hop &,
                                       const eckit::Configuration * const *,
                                       const ioda::ObsSpace &,
                                       const util::DateTime * const *,
                                       const util::DateTime * const *,
                                       const util::DateTime * const *,
                                       const util::Duration * const *);
  void ufo_timeoper_delete_f90(F90hop &);
  void ufo_timeoper_simobs_f90(const F90hop &, const F90goms &,
                               const ioda::ObsSpace &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_TIMEOPER_OBSTIMEOPER_INTERFACE_H_
