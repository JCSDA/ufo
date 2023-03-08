/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_VERTINTERP_INTERFACE_H_
#define UFO_UTILS_VERTINTERP_INTERFACE_H_

namespace ufo {

/// Interface to Fortran vertical interpolation routines

extern "C" {
void vert_interp_weights_f90(const int &nlev, const double &obl, const double *vec,
                             int &wi, double &wf);

void vert_interp_apply_f90(const int &nlev, const double *fvec,
                           double &f,
                           const int &wi, const double &wf);
void nearestneighbor_interp_index_f90(const int &nlev, const double &obl, const double *vec,
                                      int &idx);

void nearestneighbor_interp_apply_f90(const int &nlev, const double *fvec,
                                      double &f, const int &idx);
}  // extern C

}  // namespace ufo
#endif  // UFO_UTILS_VERTINTERP_INTERFACE_H_
