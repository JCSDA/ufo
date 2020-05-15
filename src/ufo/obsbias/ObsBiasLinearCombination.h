/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_OBSBIASLINEARCOMBINATION_H_
#define UFO_OBSBIAS_OBSBIASLINEARCOMBINATION_H_

#include <Eigen/Core>

#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"

#include "ufo/obsbias/ObsBiasBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class FunctionBase;
  class ObsBiasIncrement;
  class ObsDiagnostics;

/// Class to handle observation bias model using linear combination

class ObsBiasLinearCombination : public ObsBiasBase,
                                 private util::ObjectCounter<ObsBiasLinearCombination> {
 public:
  static const std::string classname() {return "ufo::ObsBiasLinearCombination";}

/// Constructor
  ObsBiasLinearCombination(const eckit::Configuration &,
                           const std::vector<std::string> &,
                           const std::vector<int> &);

/// Destructor
  ~ObsBiasLinearCombination() {}

/// I/O and diagnostics
  void read(const std::string &) override;

  void write(const eckit::Configuration &) const override;

/// diagnostics
  double norm() const override;

/// Add increments
  ObsBiasLinearCombination & operator+=(const ObsBiasIncrement &) override;
  ObsBiasLinearCombination & operator=(const ObsBias &) override;

/// Obs bias operator
  void computeObsBias(ioda::ObsVector &, const Eigen::MatrixXd &) const override;

/// save Bias Terms for QC
  void saveObsBiasTerms(ioda::ObsSpace &,
                        const std::string &,
                        const Eigen::MatrixXd &) const override;

/// Bias parameters interface
  std::size_t size() const override {return biascoeffs_.size();}
  double & operator[](const unsigned int ii) override {return biascoeffs_[ii];}

 private:
  void print(std::ostream &) const override;

  std::vector<double> biascoeffs_;
  const std::vector<std::string> prednames_;
  const std::vector<int> jobs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASLINEARCOMBINATION_H_
