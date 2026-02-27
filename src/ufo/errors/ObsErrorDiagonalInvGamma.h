/*
 * (C) Copyright 2026 Bureau of Meteorology.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORDIAGONALINVGAMMA_H_
#define UFO_ERRORS_OBSERRORDIAGONALINVGAMMA_H_

#include <memory>
#include <string>

#include "ioda/ObsVector.h"

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorParametersBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagonalInvGammaParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagonalInvGammaParameters, ObsErrorParametersBase)
 public:
  /// perturbation amplitude multiplier
  oops::Parameter<double> pert{"obs perturbations amplitude", 1.0, this};
  oops::Parameter<bool> oneMeanPerturbations{"one-mean perturbations", false, this};
  /// Relative variance parameter
  oops::RequiredParameter<double> relvar{"relative variance", this};
  /// 1-based ensemble member index.
  /// Used (and required) only if `one-mean perturbations` is set to true.
  oops::OptionalParameter<int> member{"member", this, {oops::minConstraint(1)}};
  /// Number of ensemble members.
  /// Used (and required) only if `one-mean perturbations` is set to true.
  oops::OptionalParameter<int> numberOfMembers{"number of members", this, {oops::minConstraint(1)}};
};

// -----------------------------------------------------------------------------
/// \brief Diagonal observation error covariance matrix.
class ObsErrorDiagonalInvGamma : public ObsErrorBase {
 public:
  /// The type of parameters for this class.
  typedef ObsErrorDiagonalInvGammaParameters Parameters_;

  static const std::string classname() {return "ufo::ObsErrorDiagonalInvGamma";}

  ObsErrorDiagonalInvGamma(const Parameters_ &, ioda::ObsSpace &, const eckit::mpi::Comm &);

/// Update after obs errors potentially changed
  void update(const ioda::ObsVector &) override;

/// Multiply a Departure by \f$R\f$
  void multiply(ioda::ObsVector &) const override;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ioda::ObsVector &) const override;

/// Create local R matrix
  void localize(ioda::ObsVector &) const override;

/// Return dimension of local R matrix.
  int localDim() const override;

/// Multiply a local obs vector \p zz by \f$R^{-1}\f$.
  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf &zz) const override;

/// Generate random perturbation
  void randomize(ioda::ObsVector &) const override;

/// Save obs errors
  void save(const std::string &) const override;

/// Get mean error for Jo table
  double getRMSE() const override {return stddev_.rms();}

/// Get obs errors std deviation
  std::unique_ptr<ioda::ObsVector> getObsErrors() const override;

/// Return inverseVariance
  std::unique_ptr<ioda::ObsVector> getInverseVariance() const override;

  Eigen::VectorXd local_invVarR() const {return local_inverseVariance_;}

 private:
  void print(std::ostream &) const override;
  void randomizeWithoutOneEnsembleMean(ioda::ObsVector &) const;
  void randomizeWithOneEnsembleMean(ioda::ObsVector &) const;

  ioda::ObsVector stddev_;
  ioda::ObsVector yobs_;
  ioda::ObsVector inverseVariance_;
  mutable Eigen::VectorXd local_stddev_;
  mutable Eigen::VectorXd local_inverseVariance_;
  Parameters_ params_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORDIAGONALINVGAMMA_H_
