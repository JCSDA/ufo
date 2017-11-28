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

Locations::Locations(const ObsSpace & obss,
                     const util::DateTime & t1, const util::DateTime & t2) {

  oops::Log::trace() << "Locations contructor starting " << t1 << " " << t2 << std::endl;
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  ufo_obsdb_getlocations_f90(obss.toFortran(), &p1, &p2, keyLoc_);
  oops::Log::trace() << "Locations contructor key = " << keyLoc_ << std::endl;
}

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
  os << "Locations::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace UFO

