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

#include "eckit/mpi/Comm.h"
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

  explicit Locations(const eckit::mpi::Comm &);
  explicit Locations(const ioda::ObsSpace &);
  Locations(const eckit::Configuration &, const eckit::mpi::Comm &);
  explicit Locations(const ufo::Locations &);
  ~Locations();

  Locations & operator+=(const Locations &);

  int nobs() const;
  int toFortran() const {return keyLoc_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}

 private:
  void print(std::ostream & os) const;
  F90locs keyLoc_;
  const eckit::mpi::Comm & comm_;
};

}  // namespace ufo

#endif  // UFO_LOCATIONS_H_
