/*
 * (C) Copyright 2017-2018 UCAR
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

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

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

  // Obs bias model
  void computeObsBias(ioda::ObsVector &, ObsDiagnostics &,
                      const std::vector<ioda::ObsVector> &) const;

  // Obs Bias Predictors
  std::vector<ioda::ObsVector> computePredictors(const GeoVaLs &, const ObsDiagnostics &) const;

  // Required variables
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  // Operator
  operator bool() const {return biascoeffs_.size() > 0;}

 private:
  void print(std::ostream &) const override;

  ioda::ObsSpace & odb_;

  /// bias correction coefficients (npredictors x nchannels)
  Eigen::MatrixXf biascoeffs_;

  std::vector<std::shared_ptr<PredictorBase>> predbases_;

  /// predictor names
  std::vector<std::string> prednames_;
  /// channel numbers
  std::vector<int> jobs_;

  oops::Variables geovars_;
  oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
