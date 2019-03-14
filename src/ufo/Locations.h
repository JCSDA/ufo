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

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ioda/ObsSpace.h"

#include "ufo/Locations.interface.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// Locations class to handle locations for UFO.

class Locations : public util::Printable,
                  private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "ufo::Locations";}

  Locations();
  Locations(const ioda::ObsSpace &, const util::DateTime &, const util::DateTime &);
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
