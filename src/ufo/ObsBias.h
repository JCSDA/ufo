/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_H_
#define UFO_OBSBIAS_H_

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ioda/ObsDataVector.h"

#include "ufo/ObsBiasParameters.h"
#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class ObsBiasIncrement;

/// Class to handle observation bias correction coefficients
/// \details contains information on what predictors are used for bias
///          correction application
class ObsBias : public util::Printable,
                private util::ObjectCounter<ObsBias> {
 public:
  typedef ObsBiasParameters Parameters_;

  static const std::string classname() {return "ufo::ObsBias";}

  ObsBias(ioda::ObsSpace &, const eckit::Configuration &);
  ObsBias(const ObsBias &, const bool);

  ObsBias & operator+=(const ObsBiasIncrement &);
  ObsBias & operator=(const ObsBias &);

  /// Read bias correction coefficients from file
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::size_t size() const {return biascoeffs_.size();}

  /// Return the coefficient of record \p jrec, variable \p jvar and predictor \p jpred
  ///
  /// Note: \p jpred may be the index of a static or a variable predictor.
  double operator()(size_t jrec, size_t jvar, size_t jpred) const {
    return jpred < numStaticPredictors_ ?
           1.0 : biascoeffs_[index(jrec, jvar, jpred - numStaticPredictors_)];
  }

  /// Return bias correction coefficients (for *variable* predictors)
  const Eigen::VectorXd & data() const {return biascoeffs_;}

  // Required variables
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::ObsVariables & requiredHdiagnostics() const {return hdiags_;}
  const std::vector<std::string> & requiredPredictors() const {return prednames_;}

  /// Return a reference to the vector of all (static and variable) predictors.
  const Predictors & predictors() const {return predictors_;}

  /// Return the vector of variable predictors.
  std::vector<std::shared_ptr<const PredictorBase>> variablePredictors() const;

  /// Return boolean that indicates whether different records have separate
  /// bias-correction coefficients.
  bool byRecord() const {return byRecord_;}

  /// Return the number of records that have separate bias-correction coefficients,
  /// or 1 if all obs use the same coefficients.
  const std::size_t & nrecs() const {return nrecs_;}

  /// Return the list of simulated variables.
  const oops::ObsVariables & simVars() const {return vars_;}

  /// Return the indices of variables (or channels) that don't need bias correction
  const std::vector<int> & varIndexNoBC() const {return varIndexNoBC_;}

  /// Set all variable predictors coefficients to zero (used in the test)
  void zero();

  /// \brief Initialize coefficients that have no prior value from uncorrected departures.
  ///
  /// Reproduces step [1] of WRFDA's `da_varbc_coldstart`: the coefficient of the constant
  /// predictor is set to the mode of the histogram of `y - H(x)` for each channel that is being
  /// cold-started. Does nothing unless a cold start is configured and has not yet run.
  ///
  /// \param odb      ObsSpace supplying observed values and the MPI distribution.
  /// \param hofx     Uncorrected H(x), before bias correction has been added.
  /// \param qcFlags  QC flags; only locations with flag 0 contribute.
  void coldStart(const ioda::ObsSpace & odb, const ioda::ObsVector & hofx,
                 const ioda::ObsDataVector<int> & qcFlags);

  // Operator
  operator bool() const {
    return (numStaticPredictors_ > 0 || numVariablePredictors_ > 0) && vars_.size() > 0;
  }

 private:
  void print(std::ostream &) const override;

  /// index in flattened biascoeffs_ for record \p jrec, variable \p jvar
  /// and variable predictor \p jpred
  size_t index(size_t jrec, size_t jvar, size_t jpred) const {
    return jrec * (vars_.size() * numVariablePredictors_)
             + jvar * numVariablePredictors_ + jpred;
  }

  void initPredictor(const PredictorParametersWrapper &params);


  /// bias correction coefficients (nrecords x nprimitivevariables x npredictors)
  Eigen::VectorXd biascoeffs_;

  /// bias correction predictors
  Predictors predictors_;
  /// predictor names
  std::vector<std::string> prednames_;
  /// number of static predictors (i.e. predictors with fixed coefficients all equal to 1)
  std::size_t numStaticPredictors_;
  /// number of variable predictors (i.e. predictors with variable coefficients)
  std::size_t numVariablePredictors_;

  /// store data from the input file when using by record
  /// for outputting during the write procedure
  std::vector<Eigen::ArrayXXf> inputBiasCoeffs_;
  std::vector<std::string> inputPredictors_;
  std::vector<std::string> inputRecords_;

  /// bias-correct by record?
  bool byRecord_;

  /// number of records that have separate bias-correction coefficients
  /// (nrecs_ = 1 if all obs use the same coefficients)
  std::size_t nrecs_;
  /// vector of strings of record IDs
  std::vector<std::string> recIds_;

  /// list of simulated variables
  oops::ObsVariables vars_;
  /// indices of variables (or channels) that don't need bias correction
  std::vector<int> varIndexNoBC_;

  /// Variables that need to be requested from the model (for computation of predictors)
  oops::Variables geovars_;
  /// Diagnostics that need to be requested from the obs operator (for computation of predictors)
  oops::ObsVariables hdiags_;

  /// MPI rank, used to determine whether the task should output bias coeffs to a file
  size_t rank_;

  /// MPI communicator
  const eckit::mpi::Comm & comm_;

  /// MPI communicator used in time decomposition for 4DEnVar and weak-constraint 4DVar
  const eckit::mpi::Comm & commTime_;

  /// VarBC cold-start settings, and whether the (one-shot) cold start has already run
  bool coldStartEnabled_;
  bool coldStartForce_;
  bool coldStartDone_;
  std::size_t coldStartBins_;
  double coldStartHalfWidth_;
  std::size_t coldStartMinObs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
