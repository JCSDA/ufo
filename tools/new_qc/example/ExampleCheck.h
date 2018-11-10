/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_EXAMPLECHECK_H_
#define UFO_EXAMPLECHECK_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/ExampleCheck.interface.h"

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

/// Example Check

class ExampleCheck : public util::Printable,
                     private util::ObjectCounter<ExampleCheck> {
 public:
  static const std::string classname() {return "ufo::ExampleCheck";}

  ExampleCheck(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ExampleCheck();

  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const;

 private:
  void print(std::ostream &) const;
  F90check key_;
};

}  // namespace ufo

#endif  // UFO_EXAMPLECHECK_H_
