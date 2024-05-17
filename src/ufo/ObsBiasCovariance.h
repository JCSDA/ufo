/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASCOVARIANCE_H_
#define UFO_OBSBIASCOVARIANCE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/mpi/Comm.h"

#include "oops/base/ObsVariables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/ObsBiasParameters.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ObsBias;
  class ObsBiasIncrement;
  class ObsBiasPreconditioner;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  typedef ObsBiasParameters Parameters_;

  static const std::string classname() {return "ufo::ObsBiasCovariance";}

// Constructor, destructor
  ObsBiasCovariance(ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &);
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &);
  const std::vector<std::string> predictorNames() const {return prednames_;}
  std::unique_ptr<ObsBiasPreconditioner> preconditioner() const;

 private:
  void print(std::ostream &) const {}

  ioda::ObsSpace & odb_;

// Hessian contribution from Jo bias correction terms
  std::vector<double> ht_rinv_h_;

// preconditioner
  std::vector<double> preconditioner_;

// QCed obs numbers <channel>
  std::vector<std::size_t> obs_num_;

// Minimal required QCed obs number to add contribution
  std::size_t minimal_required_obs_number_;

// Analysis error variances
  std::vector<double> analysis_variances_;

// Error variances
  Eigen::VectorXd variances_;

// Smallest variance value
  double smallest_variance_ = ObsBiasCovarianceParameters::defaultSmallestVariance();

// Largest variance value
  double largest_variance_ = ObsBiasCovarianceParameters::defaultLargestVariance();

// Largest analysis error variance
  double largest_analysis_variance_ = ObsBiasCovarianceParameters::defaultLargestAnalysisVariance();

// Step size
  double step_size_ = ObsBiasCovarianceParameters::defaultStepSize();

  std::vector<std::string> prednames_;

  /// number of records that have separate bias-correction coefficients
  /// (nrecs_ = 1 if all obs use the same coefficients)
  std::size_t nrecs_;

  /// variables for which bias correction coefficients will be updated
  oops::ObsVariables vars_;

  /// MPI rank, used to determine whether the task should output bias errors coeffs to a file
  const size_t rank_;

  /// MPI communicator used in time decomposition for 4DEnVar and weak-constraint 4DVar
  const eckit::mpi::Comm & commTime_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASCOVARIANCE_H_
