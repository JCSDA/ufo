/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLES_H_
#define UFO_VARIABLES_H_

#include <ostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "Fortran.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ufo {

// -----------------------------------------------------------------------------
/// Variables class to handle variables for the UFO

class Variables : public util::Printable,
                  private util::ObjectCounter<Variables> {
 public:
  static const std::string classname() {return "ufo::Variables";}

  explicit Variables(const oops::Variables &);
  explicit Variables(const eckit::Configuration & config);
  explicit Variables(const int keyVar);

  ~Variables();

  Variables(const Variables & other);

  F90vars & toFortran() {return keyVar_;}
  const F90vars & toFortran() const {return keyVar_;}

 private:
  void print(std::ostream & os) const;
  int keyVar_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLES_H_
