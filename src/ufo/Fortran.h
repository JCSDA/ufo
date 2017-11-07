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

/// Interface to Fortran UFO routines
/*!
 * The core of the UFO is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void ufo_var_create_f90(int & keyVars, const eckit::Configuration * const *);
  void ufo_var_clone_f90(const int & keyVars, int & keyVars_other);
  void ufo_var_info_f90(const int & keyVars, int &, int &);
  void ufo_var_delete_f90(int & keyVars);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_FORTRAN_H_
