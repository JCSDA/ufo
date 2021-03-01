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

  // Accessor to bias coefficients
  double operator[](const unsigned int ii) const {return biascoeffs_(ii);}
  double operator()(size_t ii, size_t jj) const {return biascoeffs_(ii, jj);}

  // Required variables
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  const Predictors & predictors() const {return predictors_;}
  const std::vector<int> & jobs() const {return jobs_;}

  // Operator
  operator bool() const {return biascoeffs_.size() > 0;}

 private:
  void print(std::ostream &) const override;

  /// bias correction coefficients (npredictors x nchannels)
  Eigen::MatrixXf biascoeffs_;

  /// bias correction predictors
  Predictors predictors_;

  /// predictor names
  std::vector<std::string> prednames_;
  /// channel numbers
  std::vector<int> jobs_;

  /// Variables that need to be requested from the model (for computation of predictors)
  oops::Variables geovars_;
  /// Diagnostics that need to be requested from the obs operator (for computation of predictors)
  oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
