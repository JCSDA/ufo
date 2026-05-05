/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORDIFFUSION_H_
#define UFO_ERRORS_OBSERRORDIFFUSION_H_

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
        "Gaspari-Cohn correlation lengthscale in meters", this};
  oops::RequiredParameter<int> normalizationIterations{"normalization iterations", this};
  oops::Parameter<int> InverseIterations{"Number of iterations for computing the inverse",
      100, this};
  oops::Parameter<double> InverseAccuracy{"Accuracy of inverse", 1.0e-7, this};
  oops::Parameter<bool> outputDiffusionMesh{"Output diffusion mesh", false, this};
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

  // Offset where obs locations start in the merged function space
  // (0 if no control grid, nRemainingControlPoints otherwise)
  atlas::idx_t obsOffset_ = 0;  // <-- add here

  /// Creates control grid, removes close control points and
  //  returns merged grid points (remaining control grid + obs locations)
  std::vector<atlas::PointLonLat> createControlGrid(
    const std::vector<atlas::PointLonLat>& obsnodes,
    const double removeWithin,
    const int gridSpacing);
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORDIFFUSION_H_
