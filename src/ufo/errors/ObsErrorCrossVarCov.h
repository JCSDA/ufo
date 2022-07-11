/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORCROSSVARCOV_H_
#define UFO_ERRORS_OBSERRORCROSSVARCOV_H_

#include <Eigen/Core>

#include <memory>
#include <string>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsErrorBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/ObsTraits.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for obs errors with cross-variable correlations
class ObsErrorCrossVarCovParameters : public oops::ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorCrossVarCovParameters, ObsErrorParametersBase)
 public:
  /// Input file containing correlations or covariances. If covariances are
  /// specified, they will be converted to correlations.
  oops::RequiredParameter<std::string> inputFile{"input file", this};
};

// -----------------------------------------------------------------------------
/// \brief Observation error covariance matrix with cross-variable
///        correlations.
/// \details Correlations are the same at all locations and are read
///          from the file specified in the configuration. Obs error standard
///          deviations are read from ObsSpace as ObsError group.
///          Full observation error covariance matrix is R = D^{1/2} * C * D^{1/2}
///          where D^{1/2} is a diagonal matrix with stddev_ (ObsError group)
///          on the diagonal, and C is the correlation matrix.
class ObsErrorCrossVarCov : public oops::interface::ObsErrorBase<ObsTraits> {
 public:
  /// The type of parameters passed to the constructor.
  /// This typedef is used by the ObsErrorFactory.
  typedef ObsErrorCrossVarCovParameters Parameters_;

  /// Initialize observation errors
  ObsErrorCrossVarCov(const Parameters_ &, ioda::ObsSpace &,
                      const eckit::mpi::Comm &timeComm);

  /// Update obs error standard deviations to be equal to \p stddev
  void update(const ioda::ObsVector & stddev) override;

  /// Multiply \p y by this observation error covariance
  /// Computed as R * dy = D^{1/2} * C * D^{1/2} * dy
  /// where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  ///       C - correlations
  void multiply(ioda::ObsVector & y) const override;

  /// Multiply \p y by inverse of this observation error covariance
  /// Computed as R^{-1} * dy = D^{-1/2} * C^{-1] * D^{-1/2} * dy
  /// where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  ///       C - correlations
  void inverseMultiply(ioda::ObsVector & y) const override;

  /// Generate \p y as a random perturbation
  void randomize(ioda::ObsVector & y) const override;

  /// Save obs error standard deviations under \p name group name
  void save(const std::string & name) const override;

  /// Return RMS of obs error standard deviations
  double getRMSE() const override {return stddev_.rms();}

  /// Return obs errors std deviation
  std::unique_ptr<ioda::ObsVector> getObsErrors() const override;

  /// Return inverse of obs error variance
  std::unique_ptr<ioda::ObsVector> getInverseVariance() const override;

 private:
  /// Print covariance details (for logging)
  void print(std::ostream &) const override;
  /// Observation error standard deviations
  ioda::ObsVector stddev_;
  /// Variables for which correlations are defined (same as ObsSpace::obsvariables())
  const oops::Variables vars_;
  /// Correlations between variables
  Eigen::MatrixXd varcorrelations_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORCROSSVARCOV_H_
