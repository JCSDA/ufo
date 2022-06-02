
/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_UFO_UTILS_INTERFACE_H_
#define UFO_UTILS_UFO_UTILS_INTERFACE_H_

namespace ufo {

extern "C" {
  void ufo_ops_satrad_qsplit_f90(const int &, const int &, const float *, const float *,
                                 const float *, float *, float *, float *, const bool &);
  void ufo_ops_satrad_qsatwat_f90(float *, const float *, const float *, const int &);
  void ufo_ops_qsat_f90(float *, const float *, const float *, const int &);
}  // extern C

}  // namespace ufo

#endif  //  UFO_UTILS_UFO_UTILS_INTERFACE_H_
