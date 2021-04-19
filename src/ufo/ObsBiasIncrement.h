/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <string>
#include <vector>

#include "oops/util/Printable.h"

#include "ufo/ObsBiasParameters.h"
#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ObsBias;

/// Contains increments to bias correction coefficients
class ObsBiasIncrement : public util::Printable {
 public:
  typedef ObsBiasParameters Parameters_;

  ObsBiasIncrement(const ioda::ObsSpace & odb,
                   const Parameters_ & params);
  ObsBiasIncrement(const ObsBiasIncrement &, const bool = true);

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

  /// Return the coefficient of predictor \p jpred for variable \p jvar.
  double & operator()(size_t jpred, size_t jvar) {
    return biascoeffsinc_(jvar*prednames_.size() + jpred);
  }

  /// Return bias coefficient increments
  const Eigen::VectorXd & data() const {return biascoeffsinc_;}
  Eigen::VectorXd & data() {return biascoeffsinc_;}

  // We could store coefficients in a different order. Then it would
  // be easier to extract part of the existing vector.
  std::vector<double> coefficients(size_t jpred) const {
    std::vector<double> coeffs(vars_.size());
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      coeffs[jvar] = biascoeffsinc_(jvar*prednames_.size() + jpred);
    }
    return coeffs;
  }

  // Serialize and deserialize
  std::size_t serialSize() const {return biascoeffsinc_.size();}
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, std::size_t &);

  // Operator
  operator bool() const {return biascoeffsinc_.size() > 0;}

  /// Reduce (sum) all bias correction coefficients increments
  void allSumInPlace();

 private:
  void print(std::ostream &) const;

  const ioda::ObsSpace & odb_;

  /// Bias coefficient increments
  Eigen::VectorXd biascoeffsinc_;
  std::vector<std::string> prednames_;

  /// Variables that need to be bias-corrected
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
