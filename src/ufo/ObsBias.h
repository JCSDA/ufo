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

#include "ufo/obsbias/ObsBiasBase.h"
#include "ufo/obsbias/predictors/PredictorBase.h"

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

  ObsBias(const ioda::ObsSpace &, const eckit::Configuration &);
  ObsBias(const ObsBias &, const bool);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasIncrement &);
  ObsBias & operator=(const ObsBias &);

  // I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::size_t size() const;

  // Bias parameters interface
  const double & operator[](const unsigned int ii) const {return (*biasbase_)[ii];}

  // Obs bias model
  void computeObsBias(ioda::ObsVector &, const Eigen::MatrixXd &) const;

  // Obs Bias Predictors
  Eigen::MatrixXd computePredictors(const GeoVaLs &, const ObsDiagnostics &) const;

  // Required variables
  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  // Operator
  operator bool() const {return biasbase_.get();}

 private:
  void print(std::ostream &) const;

  const ioda::ObsSpace & odb_;
  eckit::LocalConfiguration conf_;

  std::unique_ptr<ObsBiasBase> biasbase_;
  std::vector<std::shared_ptr<PredictorBase>> predbases_;
  std::vector<std::string> prednames_;
  std::vector<int> jobs_;
  oops::Variables geovars_;
  oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
