/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_
#define UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

/// Interface to Fortran UFO marine/insitutemperature routines

extern "C" {

// -----------------------------------------------------------------------------

  /// \param operatorVars
  ///   Variables to be simulated by this operator.
  /// \param operatorVarIndices
  ///   Indices of the variables from \p operatorVar in the list of all simulated
  ///   variables in the ObsSpace.
  /// \param numOperatorVarIndices
  ///   Size of the \p operatorVarIndices array (must be the same as the number of variables in
  ///   \p operatorVars).
  /// \param[out] requiredVars
  ///   GeoVaLs required for the simulation of the variables \p operatorVars.

  void ufo_insitutemperature_setup_f90(F90hop &, const eckit::Configuration &,
                                      const oops::ObsVariables &operatorVars,
                                      const int *operatorVarIndices,
                                      const int numOperatorVarIndices,
                                      oops::Variables &requiredVars);
  void ufo_insitutemperature_delete_f90(F90hop &);
  void ufo_insitutemperature_simobs_f90(const F90hop &, const F90goms &, const ioda::ObsSpace &,
                               const int &, const int &, double &);

// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_INTERFACE_H_
