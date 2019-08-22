/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_OBSBIASRADIANCEGSI_H_
#define UFO_OBSBIAS_OBSBIASRADIANCEGSI_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/obsbias/ObsBiasBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class ObsBiasIncrement;

/// Class to handle observation bias model from GSI Radiance.

class ObsBiasRadianceGSI : public ObsBiasBase,
                           private util::ObjectCounter<ObsBiasRadianceGSI> {
 public:
  static const std::string classname() {return "ufo::ObsBiasRadianceGSI";}

/// Constructor
  explicit ObsBiasRadianceGSI(const eckit::Configuration &);

/// Destructor
  virtual ~ObsBiasRadianceGSI() {}

/// I/O and diagnostics
  void read(const eckit::Configuration &) override;
  void write(const eckit::Configuration &) const override;
  double norm() const override;
  std::size_t size() const override { return bias_.size();};

/// Add increments
  ObsBiasRadianceGSI & operator+=(const ObsBiasIncrement &) override;

/// Obs bias operator
  void computeObsBias(const GeoVaLs &,
                      ioda::ObsVector &,
                      const ioda::ObsSpace &) const override;

/// Other
  const oops::Variables & variables() const override {return *varin_;}

/// Bias parameters interface
  double & operator[](const unsigned int ii) override {return bias_[ii];}
  const double & operator[](const unsigned int ii) const override {return bias_[ii];}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;
  std::string sensor_id_;  // sensor_id
  std::vector<int> channels_;  // channel

  std::vector<double> bias_;

  static const std::vector<std::string> predictors_;  // predictor names
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASRADIANCEGSI_H_
