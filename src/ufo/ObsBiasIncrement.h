/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Printable.h"

#include "ufo/predictors/PredictorBase.h"

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
// Constructor, destructor
  ObsBiasIncrement(const ioda::ObsSpace &, const eckit::Configuration &);
  ObsBiasIncrement(const ObsBiasIncrement &, const bool = true);
  ObsBiasIncrement(const ObsBiasIncrement &, const eckit::Configuration &);
  ~ObsBiasIncrement() {}

// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &);
  void zero();
  ObsBiasIncrement & operator=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator+=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator-=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator*=(const double);
  void axpy(const double, const ObsBiasIncrement &);
  double dot_product_with(const ObsBiasIncrement &) const;

// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const;

  double & operator[](const unsigned int ii) {return biascoeffsinc_[ii];}
  const double & operator[](const unsigned int ii) const {return biascoeffsinc_[ii];}

// Linear obs bias model
  void computeObsBiasTL(const GeoVaLs &,
                        const ioda::ObsDataVector<double> &,
                        ioda::ObsVector &) const;

  void computeObsBiasAD(GeoVaLs &,
                        const ioda::ObsDataVector<double> &,
                        const ioda::ObsVector &);

// Serialize and deserialize
  std::size_t serialSize() const {return biascoeffsinc_.size();}
  void serialize(std::vector<double> &) const {}
  void deserialize(const std::vector<double> &, std::size_t &) {}

// Operator
  operator bool() const {return biascoeffsinc_.size() > 0;}

 private:
  void print(std::ostream &) const;

  const ioda::ObsSpace & odb_;
  const eckit::LocalConfiguration conf_;

  std::vector<double> biascoeffsinc_;
  std::vector<std::shared_ptr<PredictorBase>> predbases_;
  std::vector<std::string> prednames_;
  std::vector<int> jobs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
