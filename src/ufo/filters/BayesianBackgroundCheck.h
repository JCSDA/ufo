/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_BAYESIANBACKGROUNDCHECK_H_
#define UFO_FILTERS_BAYESIANBACKGROUNDCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the BayesianBackgroundCheck filter.
class BayesianBackgroundCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(BayesianBackgroundCheckParameters, FilterParametersBase)

 public:
  /// Parameters related to PGE calculations
  ProbabilityOfGrossErrorParameters PGEParameters{this};
  // Probability density of observation being bad:
  oops::RequiredParameter<float> PdBad{"prob density bad obs", this,
                                       {oops::exclusiveMinConstraint(0.0f)}};
  // Initial probability of gross error:
  oops::RequiredParameter<float> PGE{"initial prob gross error", this,
                                     {oops::exclusiveMinConstraint(0.0f)}};
};

/// \brief BayesianBackgroundCheck: check observation closeness to background, accounting
///   for probability of gross error, i.e. that observation is bad.
/// \details Flag observations whose values differ from their model equivalents by too much.
/// "Too much" depends on the uncertainties of both observation and model
/// background, and the prior probability of the observation being "bad", or in
/// gross error. Can apply to scalar or vector observation types on single-level.
/// (Cannot handle profiles.) Calls function ufo::BayesianPGEUpdate.
///
/// Requires the following be specified in .yaml, under
///
/// obs filters:
/// - filter: Bayesian Background Check
///
///   * prob density bad obs: uniform probability density that observation is "bad",
///           i.e. in gross error;
///   * initial prob gross error: initial PGE value of the uniform probability
///           distribution, before adjustment depending on obs-BG closeness;
///
/// May also specify the following optional parameters (defaults in
///   ufo/utils/ProbabilityOfGrossErrorParameters.h):
///   * PGE threshold: fail if after-check PGE exceeds this value;
///   * obs minus BG threshold: fail if \f$([y-H(x)]/{\sigma})^2\f$ exceeds this;
///   * max exponent: maximum allowed value of the exponent in the 'good'
///      probability distribution;
///   * obs error multiplier: weight of observation error in combined error variance;
///   * BG error multiplier: weight of background error in combined error variance.
class BayesianBackgroundCheck : public FilterBase,
                        private util::ObjectCounter<BayesianBackgroundCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef BayesianBackgroundCheckParameters Parameters_;

  static const std::string classname() {return "ufo::BayesianBackgroundCheck";}

  BayesianBackgroundCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                  std::shared_ptr<ioda::ObsDataVector<int> > flags,
                  std::shared_ptr<ioda::ObsDataVector<float> > obserr);
  ~BayesianBackgroundCheck();

 private:
  void print(std::ostream &) const override;

  /// \brief Apply Bayesian background check filter. Return flagged=true for rejected obs.
  void applyFilter(const std::vector<bool> & apply, const Variables & filtervars,
                   std::vector<std::vector<bool>> & flagged) const override;

  /// \brief Return bayesianQC flag for observations rejected by Bayesian BG check.
  int qcFlag() const override {return QCflags::bayesianQC;}

  /// \brief Return the name of the variable containing the background error estimate of the
  /// specified filter variable.
  Variable backgrErrVariable(const Variable & filterVariable) const;

  /// \brief Reduce a vector to only the elements for which j_reduced=true.
  ///  \param[in] vector_full: full vector,
  ///  \param[in] j_reduced: indices of vector_full to be copied into...
  ///  \return vector_reduced.
  template <class T>
  std::vector<T> reduceVector(const std::vector<T> & vector_full,
                              const std::vector<size_t> & j_reduced) const;

  /// \brief Copy a reduced vector back into the correct indices of the full vector.
  ///  \param[in] vector_reduced: reduced vector to restore to...
  ///  \param[inout] vector_full: full vector,
  ///  \param[in] j_reduced: indices of vector_full to be populated.
  template <class T>
  void unreduceVector(const std::vector<T> & vector_reduced,
                      std::vector<T> & vector_full,
                      const std::vector<size_t> & j_reduced) const;
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_BAYESIANBACKGROUNDCHECK_H_
