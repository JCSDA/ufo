/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FORTRAN_H_
#define UFO_FORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace ufo {

// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Observation vector key type
typedef int F90ovec;
// Obs operator key type
typedef int F90hop;
// Observation space type
typedef int F90odb;
// Observation check key type
typedef int F90ocheck;
// Observation bias key type
typedef int F90obias;

}  // namespace ufo
#endif  // UFO_FORTRAN_H_
