/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_OBSBIASRADIANCEGSI_INTERFACE_H_
#define UFO_OBSBIAS_OBSBIASRADIANCEGSI_INTERFACE_H_

namespace ufo {

/// Interface to Fortran Obs Bias utility routines
/*!
 * The core of the utilities is coded in Fortran from GSI.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {

// -----------------------------------------------------------------------------
//  clw
// -----------------------------------------------------------------------------
  void calc_clw_f90(const int &, const float &, const float &, const int &, const int &,
                    const bool &, const bool &, const bool &, const bool &, const bool &,
                    const bool &, const bool &, const bool &, const bool &,
                    const float &, const float &, const float&, float &);
// -----------------------------------------------------------------------------

}  // extern C

}  // namespace ufo
#endif  // UFO_OBSBIAS_OBSBIASRADIANCEGSI_INTERFACE_H_
