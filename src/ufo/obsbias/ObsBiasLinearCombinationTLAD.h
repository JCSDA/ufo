/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_OBSBIASLINEARCOMBINATIONTLAD_H_
#define UFO_OBSBIAS_OBSBIASLINEARCOMBINATIONTLAD_H_

#include <memory>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/obsbias/LinearObsBiasBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class ObsBias;
  class ObsBiasIncrement;

/// Class to handle linear combination TLAD

class ObsBiasLinearCombinationTLAD : public LinearObsBiasBase,
                               private util::ObjectCounter<ObsBiasLinearCombinationTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsBiasLinearCombinationTLAD";}

/// Constructor
  ObsBiasLinearCombinationTLAD(const ioda::ObsSpace &,
                               const eckit::Configuration &,
                               const std::vector<std::string> &,
                               const std::vector<int> &);

/// Destructor
  virtual ~ObsBiasLinearCombinationTLAD() {}

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &) override;
  void zero() override;
  ObsBiasLinearCombinationTLAD & operator=(const ObsBiasIncrement &) override;
  ObsBiasLinearCombinationTLAD & operator+=(const ObsBiasIncrement &) override;
  ObsBiasLinearCombinationTLAD & operator-=(const ObsBiasIncrement &) override;
  ObsBiasLinearCombinationTLAD & operator*=(const double) override;
  void axpy(const double, const ObsBiasIncrement &) override;
  double dot_product_with(const ObsBiasIncrement &) const override;

/// I/O and diagnostics
  void read(const eckit::Configuration &) const override;
  void write(const eckit::Configuration &) const override;
  double norm() const override;

/// Linear obs bias operator
  void computeObsBiasTL(const GeoVaLs &,
                        const Eigen::MatrixXd &,
                        ioda::ObsVector &) const override;

  void computeObsBiasAD(GeoVaLs &,
                        const Eigen::MatrixXd &,
                        const ioda::ObsVector &) override;

/// Bias parameters interface
  double & operator[](const unsigned int ii) override {return biascoeffsinc_[ii];}
  const double & operator[](const unsigned int ii) const override {return biascoeffsinc_[ii];}
  const std::vector<double>::iterator begin() {return biascoeffsinc_.begin();}
  const std::vector<double>::iterator end() {return biascoeffsinc_.end();}

  const ioda::ObsSpace & obsspace() const override {return odb_;}

 private:
  void print(std::ostream &) const override;

  const ioda::ObsSpace & odb_;

  std::vector<double> biascoeffsinc_;
  const std::vector<std::string> prednames_;
  const std::vector<int> jobs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASLINEARCOMBINATIONTLAD_H_
