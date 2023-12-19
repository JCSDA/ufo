/*
 * (C) Copyright 2020 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_INTERFACE_H_
#define UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_INTERFACE_H_

#include "oops/base/Variables.h"
#include "ufo/Fortran.h"

namespace ufo {

extern "C" {

// call Fortran code to get ESG grid information

  void lam_domaincheck_esg_f90(const float &, const float &, const float &, const float &,
                               const float &, const int &, const int &,
                               const float &, const float &, const int &, const float &,
                               const float &, int &);

  void lam_domaincheck_circle_f90(const float &, const float &, const float &,
                                  const float &, const float &, int &);

}  // extern C

}  // namespace ufo
#endif  // UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_INTERFACE_H_
