/*
 * (C) Copyright 2017-2024 UCAR
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASINCREMENT_H_
#define UFO_OBSBIASINCREMENT_H_

#include <Eigen/Core>

#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/base/ObsVariables.h"
#include "oops/util/Printable.h"

#include "ufo/ObsBiasParameters.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ObsBias;

/// Contains increments to bias correction coefficients
class ObsBiasIncrement : public util::Printable {
 public:
  typedef ObsBiasParameters Parameters_;

  static const std::string classname() {return "ufo::ObsBiasIncrement";}

  ObsBiasIncrement(const ioda::ObsSpace & odb, const eckit::Configuration &);

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
  void read(const eckit::Configuration &);  // For tests only, due to lack of appropriate use case;
                                            // hence it checks for YAML key "test input file"
  void write(const eckit::Configuration &) const;
  double norm() const;

  /// Return bias coefficient increments
  const Eigen::VectorXd & data() const {return biascoeffsinc_;}
  Eigen::VectorXd & data() {return biascoeffsinc_;}

  /// Return bias coefficient increments for one predictor
  std::vector<double> coefficients(size_t) const;
  /// Increment bias coeffiecient increments for one predictor
  void updateCoeff(size_t, const std::vector<double> &);

  // Serialize and deserialize
  std::size_t serialSize() const {return biascoeffsinc_.size();}
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, std::size_t &);

  // Operator
  operator bool() const {return biascoeffsinc_.size() > 0;}

 private:
  void print(std::ostream &) const;

  /// index in flattened biascoeffsinc_ for record \p jrec, variable \p jvar
  /// and variable predictor \p jpred
  size_t index(size_t jrec, size_t jvar, size_t jpred) const {
    return jrec * (vars_.size() * prednames_.size())
             + jvar * prednames_.size() + jpred;
  }

  /// Bias coefficient increments
  Eigen::VectorXd biascoeffsinc_;
  std::vector<std::string> prednames_;

  /// Bias-correct by record?
  bool byRecord_;

  /// Number of records that have separate bias-correction coefficients
  /// (nrecs_ = 1 if all obs use the same coefficients)
  std::size_t nrecs_;

  /// Vector of strings of record IDs
  std::vector<std::string> recIds_;

  /// List of simulated variables
  oops::ObsVariables vars_;

  /// Name of the output file of the bias coeff increments
  std::string outputFile_;

  /// MPI rank, used to determine whether the task should output bias coeff increments to a file
  size_t rank_;

  /// MPI communicator used in time decomposition for 4DEnVar and weak-constraint 4DVar
  const eckit::mpi::Comm & commTime_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASINCREMENT_H_
