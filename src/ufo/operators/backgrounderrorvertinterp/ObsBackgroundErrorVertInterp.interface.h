/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_INTERFACE_H_
#define UFO_OPERATORS_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_INTERFACE_H_

#include "ufo/Fortran.h"

namespace oops {
class Variables;
}  // namespace oops

namespace ioda {
class ObsSpace;
}  // namespace ioda

namespace ufo {

extern "C" {

  void ufo_backgrounderrorvertinterp_fillobsdiags_f90(const int &len_obs_vcoord,
                                                      const char *obs_vcoord,
                                                      const int &len_obs_vgroup,
                                                      const char *obs_vgroup,
                                                      const int &len_vcoord,
                                                      const char *vcoord,
                                                      const F90goms &geovals,
                                                      const ioda::ObsSpace &obsspace,
                                                      const int &nlocs,
                                                      const oops::Variables &obsvars,
                                                      const F90goms &obsdiags);

}  // extern C

}  // namespace ufo

#endif  // UFO_OPERATORS_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_INTERFACE_H_
