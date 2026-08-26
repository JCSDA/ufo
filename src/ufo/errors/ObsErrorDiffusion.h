/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORDIFFUSION_H_
#define UFO_ERRORS_OBSERRORDIFFUSION_H_

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <vector>

#include "atlas/functionspace.h"

#include "ioda/ObsVector.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorParametersBase.h"
#include "ufo/ObsTraits.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for correlated obs errors
class ObsErrorDiffusionParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiffusionParameters, ObsErrorParametersBase)

 public:
  class ControlGrid : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ControlGrid, oops::Parameters)
   public:
    oops::RequiredParameter<int> gridSpacing{"grid spacing in degrees", this};
    oops::RequiredParameter<double> removeWithin{"remove within meters", this};
  };

  oops::RequiredParameter<std::vector<std::string>> var{"correlation variable names",
        "Group/Name obs variable that correlations should be computed for (note: "
        "this variable should be the same variable as obs space is grouped on", this};
  oops::RequiredParameter<double> lscale{"correlation lengthscale",
        "Correlation lengthscale in meters, interpreted by default as the "
        "Gaspari-Cohn cutoff radius. See `as gaussian`.", this};
  oops::Parameter<bool> asGaussian{"as gaussian",
        "If false (default), interpret `correlation lengthscale` as a Gaspari-Cohn "
        "cutoff radius. If true, pass the value through as a Daley length. Mirrors the "
        "SABER diffusion block's flag.", false, this};
  oops::RequiredParameter<int> normalizationIterations{"normalization iterations", this};
  oops::Parameter<int> InverseIterations{"Number of iterations for computing the inverse",
      100, this};
  oops::Parameter<double> InverseAccuracy{"Accuracy of inverse", 1.0e-7, this};
  oops::Parameter<double> diagonalLoading{"diagonal loading",
        "Nugget fraction α in [0, 1]: blends a diagonal (uncorrelated) component "
        "into R. α=0.0 is pure correlated, α=1.0 pure diagonal. Represents the "
        "uncorrelated variance fraction; α>0.0 also speeds GMRESR R^-1 convergence "
        "at wide L_R/Δx.", 0.0, this,
        {oops::minConstraint(0.0), oops::maxConstraint(1.0)}};
  oops::Parameter<bool> allowAnyDistribution{"allow any distribution",
        "Allow a non-Atlas distribution at more than one PE. Only set this if it is known to "
        "be geographically contiguous, as the diffusion halo requires.", false, this};
  oops::Parameter<bool> outputDiffusionMesh{"output diffusion mesh", false, this};
  // Add parameter for control grid
  oops::OptionalParameter<ControlGrid> controlGrid{"control grid", this};
};


// -----------------------------------------------------------------------------
/// \brief Correlated observation error covariance matrix.
class ObsErrorDiffusion : public ObsErrorBase {
 public:
  /// The type of parameters passed to the constructor.
  /// This typedef is used by the ObsErrorFactory.
  typedef ObsErrorDiffusionParameters Parameters_;
  static const std::string classname() {return "ufo::ObsErrorDiffusion";}

  /// Initialize observation errors
  ObsErrorDiffusion(const Parameters_ &, ioda::ObsSpace &,
                    const eckit::mpi::Comm &);

  /// Update obs error standard deviations to be equal to \p stddev
  void update(const ioda::ObsVector & stddev) override;

  /// Multiply \p y by this observation error covariance
  /// Computed as R * dy = D^{1/2} * C * D^{1/2} * dy
  /// where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  ///       C - correlations
  void multiply(ioda::ObsVector & y) const override;

  /// Wrapper method needed for GMRESR takes an obs vector
  ///  \p y, multiplies it by the ObsError, and returns the
  /// result as \p y_o
  void multiply(const ioda::ObsVector & y, ioda::ObsVector & y_o) const;

  /// Multiply \p y by inverse of this observation error covariance
  /// Computed as R^{-1} * dy = D^{-1/2} * C^{-1] * D^{-1/2} * dy
  /// where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  ///       C - correlations
  void inverseMultiply(ioda::ObsVector & y) const override;

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

  /// Create local R matrix
  void localize(ioda::ObsVector &) const override;

  /// Return dimension of local R matrix.
  int localDim() const override;

  /// Multiply a local obs vector \p zz by \f$R^{-1}\f$.
  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf &zz) const override;

  Eigen::VectorXd local_invVarR() const;

 private:
  /// Print covariance details
  void print(std::ostream &) const override;
  std::shared_ptr<oops::GeometryData> geom_;
  std::shared_ptr<oops::Diffusion::DerivedGeom> diffusionGeom_;

  /// Observation error standard deviations
  Parameters_ params_;
  ioda::ObsVector stddev_;
  ioda::ObsVector inverseVariance_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<oops::Diffusion> diffusion_;
  atlas::Field hzNorm_;

  // Mesh + normalization are built once on the first update() and reused on later
  // calls (see the reuse guard in update() for why and the caveats).
  bool meshAndNormBuilt_ = false;
  int builtObsCount_ = 0;  // this rank's obs count at build time; reuse tripwire

  // Index where obs start in the global grid numbering; control points (if any)
  // fill [0, obsOffset_).
  atlas::idx_t obsOffset_ = 0;

  // Per local ObsVector location: node index in this PE's mesh for each valid obs,
  // or -1 for QC-failed locations. Used by multiply()/randomize() to copy values
  // to/from the diffusion field.
  std::vector<atlas::idx_t> localObs2node_;

  /// Build the global control grid from obs locations, dropping points within
  /// \p removeWithin of any obs; each kept point inherits its nearest obs's owner.
  void createControlGrid(
    const std::vector<atlas::PointLonLat>& obsnodes,
    const std::vector<int>& obsOwner,
    const double removeWithin,
    const int gridSpacing,
    std::vector<atlas::PointLonLat>& ctrlPoints,
    std::vector<int>& ctrlOwner) const;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORDIFFUSION_H_
