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

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/ObsBiasParameters.h"
#include "ufo/predictors/PredictorBase.h"

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVals;
  class ObsBiasIncrement;
  class ObsDiagnostics;

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

  /// Return the coefficient of predictor \p jpred for variable \p jvar.
  ///
  /// Note: \p jpred may be the index of a static or a variable predictor.
  double operator()(size_t jpred, size_t jvar) const {
    return jpred < numStaticPredictors_ ? 1.0 : biascoeffs_(jpred - numStaticPredictors_, jvar);
  }

  /// Return the ii'th element of the flattened array of coefficients of *variable* predictors.
  double operator[](const unsigned int ii) const {return biascoeffs_(ii);}

  // Required variables
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}
  const std::vector<std::string> & requiredPredictors() const {return prednames_;}

  /// Return a reference to the vector of all (static and variable) predictors.
  const Predictors & predictors() const {return predictors_;}

  /// Return the vector of variable predictors.
  std::vector<std::shared_ptr<const PredictorBase>> variablePredictors() const;

  /// Return the list of bias-corrected variables.
  const oops::Variables & correctedVars() const {return vars_;}

  // Operator
  operator bool() const {
    return (numStaticPredictors_ > 0 || numVariablePredictors_ > 0) && vars_.size() > 0;
  }

 private:
  void print(std::ostream &) const override;

  void initPredictor(const eckit::Configuration &predictorConf);

  /// bias correction coefficients (npredictors x nprimitivevariables)
  Eigen::MatrixXf biascoeffs_;

  /// bias correction predictors
  Predictors predictors_;
  /// predictor names
  std::vector<std::string> prednames_;
  /// number of static predictors (i.e. predictors with fixed coefficients all equal to 1)
  std::size_t numStaticPredictors_;
  /// number of variable predictors (i.e. predictors with variable coefficients)
  std::size_t numVariablePredictors_;

  /// corrected variables names (for now has to be same as "simulated variables")
  oops::Variables vars_;

  /// Variables that need to be requested from the model (for computation of predictors)
  oops::Variables geovars_;
  /// Diagnostics that need to be requested from the obs operator (for computation of predictors)
  oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
