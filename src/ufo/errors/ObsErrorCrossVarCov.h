/*
 * (C) Copyright 2021 UCAR.
 * (C) Crown copyright 2023 Met Office
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

#include "oops/base/ObsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorParametersBase.h"
#include "ufo/errors/ObsErrorReconditioner.h"
#include "ufo/ObsTraits.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for obs errors with cross-variable correlations
class ObsErrorCrossVarCovParameters : public ObsErrorParametersBase {
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
///          on the diagonal, and C is the correlation matrix. The cross-variable
///          R matrices at each location can be reconditioned to reduce their condition
///          number, in order to speed up convergence of minimisation.
class ObsErrorCrossVarCov : public ObsErrorBase {
 public:
  /// The type of parameters for this class.
  typedef ObsErrorCrossVarCovParameters Parameters_;

  /// Initialize observation errors
  ObsErrorCrossVarCov(const Parameters_ &, ioda::ObsSpace &,
                      const eckit::mpi::Comm &);

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

  /// Create local R matrix
  void localize(ioda::ObsVector & locvector) const override;

  /// Return dimension of local R matrix.
  int localDim() const override;

  /// Multiply a local obs vector \p zz by \f$R^{-1}\f$.
  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf &zz) const override;

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

  Eigen::VectorXd local_invVarR() const;

 private:
  /// Print covariance details (for logging)
  void print(std::ostream &) const override;
  /// Recondition the R matrix - called by update
  void recondition(const ioda::ObsVector & mask);
  /// Configuration as a Parameters_
  Parameters_ params_;
  /// Observation error standard deviations
  ioda::ObsVector stddev_;
  /// Variables for which correlations are defined (same as ObsSpace::obsvariables())
  const oops::ObsVariables vars_;
  /// Correlations between variables
  Eigen::MatrixXd varcorrelations_;
  /// Create reconditioner
  std::unique_ptr<ObsErrorReconditioner> reconditioner_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORCROSSVARCOV_H_
