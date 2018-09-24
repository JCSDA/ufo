/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSCHECK_H_
#define UFO_OBSCHECK_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "FortranObsCheck.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// ObsCheck: check observation for quality

class ObsCheck : public util::Printable,
                private util::ObjectCounter<ObsCheck> {
 public:
  static const std::string classname() {return "ufo::ObsCheck";}

  ObsCheck(const ioda::ObsSpace &, const oops::Variables &,
          const util::DateTime &, const util::DateTime &);
  explicit ObsCheck(const ioda::ObsSpace &);
  explicit ObsCheck(const eckit::Configuration &);

  ~ObsCheck();

  void postFilter(const GeoVaLs &, const ioda::ObsVector &, const ioda::ObsSpace &) const;
  void priorFilter(const ioda::ObsSpace &) const;

  int & toFortran() {return keyObsCheck_;}
  const int & toFortran() const {return keyObsCheck_;}

 private:
  void print(std::ostream &) const;
  F90ocheck keyObsCheck_;
};

}  // namespace ufo

#endif  // UFO_OBSCHECK_H_
