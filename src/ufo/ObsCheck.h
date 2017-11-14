/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSCHECK_H_
#define UFO_OBSCHECK_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class ObsSpace;
  class Variables;
  class GeoVaLs;
  class ObsVector;

/// ObsCheck: check observation for quality

class ObsCheck : public util::Printable,
                private util::ObjectCounter<ObsCheck> {
 public:
  static const std::string classname() {return "ufo::ObsCheck";}

  ObsCheck(const ObsSpace &, const Variables &,
          const util::DateTime &, const util::DateTime &);
  ObsCheck(const ObsSpace &);
  ObsCheck(const eckit::Configuration &);

  ~ObsCheck();

  void postFilter(const GeoVaLs &, const ObsVector &, const ObsSpace &) const;

  int & toFortran() {return keyObsCheck_;}
  const int & toFortran() const {return keyObsCheck_;}

 private:
  void print(std::ostream &) const;
  F90ocheck keyObsCheck_;
};

}  // namespace ufo

#endif  // UFO_OBSCHECK_H_
