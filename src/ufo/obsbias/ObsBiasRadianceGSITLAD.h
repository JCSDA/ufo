/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_OBSBIASRADIANCEGSITLAD_H_
#define UFO_OBSBIAS_OBSBIASRADIANCEGSITLAD_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/obsbias/LinearObsBiasBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class ObsBias;
  class ObsBiasIncrement;

/// Class to handle linear observation bias model from GSI Radiance.

class ObsBiasRadianceGSITLAD : public LinearObsBiasBase,
                               private util::ObjectCounter<ObsBiasRadianceGSITLAD> {
 public:
  static const std::string classname() {return "ufo::ObsBiasRadianceGSITLAD";}

/// Constructor
  explicit ObsBiasRadianceGSITLAD(const eckit::Configuration &);

/// Destructor
  virtual ~ObsBiasRadianceGSITLAD() {}

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &) override;
  void zero() override;
  void random() override;
  ObsBiasRadianceGSITLAD & operator=(const ObsBiasIncrement &) override;
  ObsBiasRadianceGSITLAD & operator+=(const ObsBiasIncrement &) override;
  ObsBiasRadianceGSITLAD & operator-=(const ObsBiasIncrement &) override;
  ObsBiasRadianceGSITLAD & operator*=(const double) override;
  void axpy(const double, const ObsBiasIncrement &) override;
  double dot_product_with(const ObsBiasIncrement &) const override;

/// I/O and diagnostics
  void read(const eckit::Configuration &) override;
  void write(const eckit::Configuration &) const override;
  double norm() const override;
  std::size_t size() const override { return biascoeffsinc_.size();};

/// Linear obs bias operator
  void computeObsBiasTL(const GeoVaLs &,
                        ioda::ObsVector &,
                        const ioda::ObsSpace &) const override;

  void computeObsBiasAD(GeoVaLs &,
                        const ioda::ObsVector &,
                        const ioda::ObsSpace &) override;

/// Other
  const oops::Variables & variables() const override {return *varin_;}

/// Bias parameters interface
  double & operator[](const unsigned int ii) override {return biascoeffsinc_[ii];}
  const double & operator[](const unsigned int ii) const override {return biascoeffsinc_[ii];}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;
  std::string sensor_id_;  // sensor_id
  std::vector<int> channels_;  // channel

  std::vector<double> biascoeffsinc_;
  std::vector<std::string> predictors_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASRADIANCEGSITLAD_H_
