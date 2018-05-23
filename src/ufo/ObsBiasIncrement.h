/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <iostream>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class ObsBias;

// -----------------------------------------------------------------------------

class ObsBiasIncrement : public util::Printable {
 public:
/// Constructor, destructor
  explicit ObsBiasIncrement() {}
  explicit ObsBiasIncrement(const eckit::Configuration &) {}
  ObsBiasIncrement(const ObsBiasIncrement &, const bool copy = true) {}
  ObsBiasIncrement(const ObsBiasIncrement &, const eckit::Configuration &) {}
  ~ObsBiasIncrement() {}

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &) {}
  void zero() {}
  ObsBiasIncrement & operator=(const ObsBiasIncrement &) {}
  ObsBiasIncrement & operator+=(const ObsBiasIncrement &) {}
  ObsBiasIncrement & operator-=(const ObsBiasIncrement &) {}
  ObsBiasIncrement & operator*=(const double) {}
  void axpy(const double, const ObsBiasIncrement &) {}
  double dot_product_with(const ObsBiasIncrement &) const {
    oops::Log::trace() << "ufo::ObsBiasIncrement dot product" << std::endl;
    std::string ishouldnotexist="If you delete me, Jb Obs Aux will be set to huge number (singularity, building in debug mode)";
    std::cout << "Ouch: " << ishouldnotexist << std::endl;
  }

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream &) const {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
