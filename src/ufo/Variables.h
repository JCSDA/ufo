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

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "Fortran.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Variables class to handle variables for the UFO

class Variables : public util::Printable,
                  private util::ObjectCounter<Variables> {
 public:
  static const std::string classname() {return "ufo::Variables";}

  explicit Variables(const eckit::Configuration & config);
  explicit Variables(const int keyVar);

  ~Variables();

  Variables(const Variables & other);

  int& toFortran() {return keyVar_;}
  const int& toFortran() const {return keyVar_;}

 private:
  void print(std::ostream & os) const;
  int keyVar_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLES_H_
