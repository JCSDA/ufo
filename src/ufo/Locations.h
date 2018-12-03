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

#include "Fortran.h"
#include "FortranLocations.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// Locations class to handle locations for UFO.

class Locations : public util::Printable,
                  private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "ufo::Locations";}

  explicit Locations(const F90locs key): keyLoc_(key) {}
  explicit Locations(const eckit::Configuration &);
  ~Locations();

  int nobs() const;
  int toFortran() const {return keyLoc_;}

 private:
  void print(std::ostream & os) const;
  F90locs keyLoc_;
};

}  // namespace ufo

#endif  // UFO_LOCATIONS_H_
