/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "Locations.h"
#include "Fortran.h"

namespace ufo {

// -----------------------------------------------------------------------------

Locations::~Locations() {
  ufo_locs_delete_f90(keyLoc_);
}

// -----------------------------------------------------------------------------

int Locations::nobs() const {
  int nobs;
  ufo_locs_nobs_f90(keyLoc_, nobs);
  return nobs;
}

// -----------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  int nobs;
  ufo_locs_nobs_f90(keyLoc_, nobs);
  os << "Locations: " << nobs << " locations";
}

// -----------------------------------------------------------------------------

}  // namespace UFO

