/*
 * (C) Copyright 2017-2021 UCAR
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

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/ObsBiasParameters.h"
#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
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

  ObsBias(ioda::ObsSpace &, const Parameters_ &);
  ObsBias(const ObsBias &, const bool);

  ObsBias & operator+=(const ObsBiasIncrement &);
  ObsBias & operator=(const ObsBias &);

  /// Read bias correction coefficients from file
  void read(const Parameters_ &);
  void write(const Parameters_ &) const;
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
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}
  const std::vector<std::string> & requiredPredictors() const {return prednames_;}

  /// Return a reference to the vector of all (static and variable) predictors.
  const Predictors & predictors() const {return predictors_;}

  /// Return the vector of variable predictors.
  std::vector<std::shared_ptr<const PredictorBase>> variablePredictors() const;

  /// Return the number of records that are bias-corrected independently from each other,
  /// or 1 if all obs are bias-corrected together
  const std::size_t & nrecs() const {return nrecs_;}

  /// Return the list of simulated variables.
  const oops::Variables & simVars() const {return vars_;}

  /// Return the indices of variables (or channels) that don't need bias correction
  const std::vector<int> & varIndexNoBC() const {return varIndexNoBC_;}

  /// Set all variable predictors coefficients to zero (used in the test)
  void zero();

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

  /// number of records that are bias-corrected independently from each other
  /// (nrecs_ = 1 if all obs are bias-corrected together)
  std::size_t nrecs_;

  /// list of simulated variables
  oops::Variables vars_;
  /// indices of variables (or channels) that don't need bias correction
  std::vector<int> varIndexNoBC_;

  /// Variables that need to be requested from the model (for computation of predictors)
  oops::Variables geovars_;
  /// Diagnostics that need to be requested from the obs operator (for computation of predictors)
  oops::Variables hdiags_;

  /// MPI rank, used to determine whether the task should output bias coeffs to a file
  size_t rank_;

  /// MPI communicator used in time decomposition for 4DEnVar and weak-constraint 4DVar
  const eckit::mpi::Comm & commTime_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
