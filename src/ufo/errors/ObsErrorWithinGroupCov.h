/*
 * (C) Copyright 2023- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsVector.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorParametersBase.h"
#include "ufo/errors/ObsErrorReconditioner.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

enum class DistanceFunctions {
      LINEAR,
      HAVERSINE
};

struct DistanceFunctionsParameterTraitsHelper {
  typedef DistanceFunctions EnumType;
  static constexpr char enumTypeName[] = "DistanceFunctions";
  static constexpr util::NamedEnumerator<DistanceFunctions> namedValues[] = {
    { DistanceFunctions::LINEAR, "linear" },
    { DistanceFunctions::HAVERSINE, "haversine" }
  };
};

enum class CorrelationFunctions {
      GC99,
      MARKOV,
      GAUSSIAN
};

struct CorrelationFunctionsParameterTraitsHelper {
  typedef CorrelationFunctions EnumType;
  static constexpr char enumTypeName[] = "CorrelationFunctions";
  static constexpr util::NamedEnumerator<CorrelationFunctions> namedValues[] = {
    { CorrelationFunctions::GC99, "gc99" },
    { CorrelationFunctions::MARKOV, "markov" },
    { CorrelationFunctions::GAUSSIAN, "gaussian"}
  };
};

}  // namespace ufo

namespace oops {

// Add the specialization for DistanceFunctions:
template <>
struct ParameterTraits<ufo::DistanceFunctions> :
    public EnumParameterTraits<ufo::DistanceFunctionsParameterTraitsHelper>
{};

// Add the specialization for CorrelationFunctions:
template <>
struct ParameterTraits<ufo::CorrelationFunctions> :
    public EnumParameterTraits<ufo::CorrelationFunctionsParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Parameters for obs errors with correlations between obs in one group
///        set by obs space.obsgrouping. The parameters are listed in alphabetical
///        order based on the parameter name.
class ObsErrorWithinGroupCovParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorWithinGroupCovParameters, ObsErrorParametersBase)
 public:
  oops::Parameter<CorrelationFunctions> correlationFunction{"correlation function",
        "Correlation function to use for correlation computations. "
        "Currently only 'gc99', 'markov', and 'gaussian' are supported",
        CorrelationFunctions::GC99, this};
  oops::RequiredParameter<double> lscale{"correlation lengthscale",
        "Correlation lengthscale used with the correlation functions", this};
  oops::OptionalParameter<std::vector<std::string>> var{"correlation variable names",
        "Group/Name obs variable that correlations should be computed for (note: "
        "this variable should be the same variable as obs space is grouped on", this};
  oops::Parameter<DistanceFunctions> distanceFunction{"distance function",
        "Distance function to use for correlation computations. "
        "Currently only 'linear' and 'haversine' are supported",
        DistanceFunctions::LINEAR, this};
  oops::Parameter<double> markovLengthscaleFactor{"lengthscale factor for markov correlation limit",
        "The lengthscale factor is multiplied by the lengthscale to provide the limit in "
        "which correlations are evaluated. Beyond this distance correlation values are set to "
        "zero.", 1.0, this};
  oops::Parameter<bool> applyBasicReconditioning{"apply basic reconditioning",
        "Apply a simple ridge regression reconditiong to the gaussian correlation profile."
        "This is the same reconditioning that is applied to the markov correlation method",
        false, this};
};

// -----------------------------------------------------------------------------
/// \brief Observation error covariance matrix with correlations between
///        obs in one group (set by obs space.obsgrouping)
/// \details Correlations are computed as GC(x1-x2) where x1, x2 are values of
///          Group/Name variable `var`. Only correlations between locations in
///          the same record are considered, locations in different records are
///          considered uncorrelated. Different variables are considered uncorrelated.
///          Obs error standard deviations are read from ObsSpace as ObsError group.
///          Full observation error covariance matrix is R = D^{1/2} * C * D^{1/2}
///          where D^{1/2} is a diagonal matrix with stddev_ (ObsError group)
///          on the diagonal, and C is the correlation matrix.
class ObsErrorWithinGroupCov : public ObsErrorBase {
 public:
  /// The type of parameters passed to the constructor.
  /// This typedef is used by the ObsErrorFactory.
  typedef ObsErrorWithinGroupCovParameters Parameters_;

  /// Initialize observation errors
  ObsErrorWithinGroupCov(const Parameters_ &, ioda::ObsSpace &,
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

  /// Save the correlations for record \p irec in the file with
  /// \p filename for diagnostics, as well as random vector \p randomVec
  /// and \p randomVec multiplied by the correlation matrix
  void saveCorrelations(const std::string & filename, size_t irec,
                        ioda::ObsVector & randomVec) const;

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
  /// Multiply only by correlations (diagnostics)
  void multiplyCorrelations(ioda::ObsVector & y) const;
  /// Recondition the R matrix - called by update
  void recondition(const ioda::ObsVector & mask);
  /// ObsError configuration as oops::Paramters
  Parameters_ params_;
  /// ObsSpace reference
  const ioda::ObsSpace & obspace_;
  /// Coordinate used for correlation computations
  ioda::ObsDataVector<float> coord_;
  /// Observation error standard deviations
  ioda::ObsVector stddev_;
  /// Each element holds a lower-triangle of the correlation matrix for all locations
  /// within one group
  std::vector<Eigen::MatrixXd> correlations_;
  /// Create reconditioner
  std::unique_ptr<ObsErrorReconditioner> reconditioner_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
