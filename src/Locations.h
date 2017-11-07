/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_LOCATIONS_H_
#define UFO_LOCATIONS_H_

#include <ostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "Fortran.h"
#include "ObsSpace.h"

namespace ufo {

/// Locations class to handle locations for the UFO

class Locations : public util::Printable,
                  private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "ufo::Locations";}

  Locations(const ObsSpace & ot,
            const util::DateTime & t1, const util::DateTime & t2) {
    keyLoc_ = ot.locations(t1, t2);
  }

  ~Locations() {qg_loc_delete_f90(keyLoc_);}

  int nobs() const {
    int nobs;
    qg_loc_nobs_f90(keyLoc_, nobs);
    return nobs;
  }

  int toFortran() const {return keyLoc_;}
 private:
  void print(std::ostream & os) const {
    os << "ufo::Locations::print not implemented";
  }
  F90locs keyLoc_;
};

}  // namespace ufo

#endif  // UFO_LOCATIONS_H_
