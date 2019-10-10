/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <iostream>
#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Printable.h"

#include "ufo/obsbias/LinearObsBiasBase.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
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
  std::size_t size() const;

  double & operator[](const unsigned int ii) {return (*biasbase_)[ii];}
  const double & operator[](const unsigned int ii) const {return (*biasbase_)[ii];}

/// Linear obs bias model
  void computeObsBiasTL(const GeoVaLs &,
                        ioda::ObsVector &,
                        const ioda::ObsSpace &) const;

  void computeObsBiasAD(GeoVaLs &,
                        const ioda::ObsVector &,
                        const ioda::ObsSpace &);

/// Serialize and deserialize
  std::size_t serialSize() const {return 0;}
  void serialize(std::vector<double> &) const {}
  void deserialize(const std::vector<double> &, std::size_t &) {}

/// Other
  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<LinearObsBiasBase> biasbase_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
