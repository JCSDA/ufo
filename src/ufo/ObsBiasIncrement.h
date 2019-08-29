/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <iostream>
#include <vector>

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
  ObsBiasIncrement() {}
  explicit ObsBiasIncrement(const eckit::Configuration &);
  ObsBiasIncrement(const ObsBiasIncrement &, const bool = true);
  ObsBiasIncrement(const ObsBiasIncrement &, const eckit::Configuration &);
  ~ObsBiasIncrement() {}

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &);
  void zero();
  void random();
  ObsBiasIncrement & operator=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator+=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator-=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator*=(const double);
  void axpy(const double, const ObsBiasIncrement &);
  double dot_product_with(const ObsBiasIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const;
  std::size_t size() const {return biasinc_.size();}

  double & operator[](const unsigned int ii) {return biasinc_[ii];}
  const double & operator[](const unsigned int ii) const {return biasinc_[ii];}

/// Serialize and deserialize
  size_t serialSize() const {return 0;}
  void serialize(std::vector<double> &) const {}
  void deserialize(const std::vector<double> &, size_t &) {}

 private:
  void print(std::ostream &) const {}
  std::vector<double> biasinc_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
