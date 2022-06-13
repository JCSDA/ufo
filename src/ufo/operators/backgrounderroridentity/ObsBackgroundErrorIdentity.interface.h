/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_INTERFACE_H_
#define UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_INTERFACE_H_

#include "ufo/Fortran.h"

namespace oops {
class Variables;
}  // namespace oops

namespace ioda {
class ObsSpace;
}  // namespace ioda

namespace ufo {

extern "C" {

  void ufo_backgrounderroridentity_fillobsdiags_f90(const F90goms &geovals,
                                                    const int &nlocs,
                                                    const oops::Variables &obsvars,
                                                    const F90goms &obsdiags);

}  // extern C

}  // namespace ufo

#endif  // UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_INTERFACE_H_
