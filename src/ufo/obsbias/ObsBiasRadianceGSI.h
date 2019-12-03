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
  class ObsDiagnostics;

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
  std::size_t size() const override { return biascoeffs_.size();}

/// Add increments
  ObsBiasRadianceGSI & operator+=(const ObsBiasIncrement &) override;
  ObsBiasRadianceGSI & operator=(const ObsBias &) override;

/// Obs bias operator
  void computeObsBias(const GeoVaLs &,
                      ioda::ObsVector &,
                      const ioda::ObsSpace &,
                      const ObsDiagnostics &) const override;

/// Obs bias predictor
  void computeObsBiasPredictors(const GeoVaLs &,
                                const ioda::ObsSpace &,
                                const ObsDiagnostics &,
                                std::vector<float> &) const;

/// Other
  const oops::Variables & requiredGeoVaLs() const override {return *geovars_;}
  const oops::Variables & requiredHdiagnostics() const override {return *hdiags_;}

/// Bias parameters interface
  double & operator[](const unsigned int ii) override {return biascoeffs_[ii];}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> geovars_;
  std::unique_ptr<const oops::Variables> hdiags_;
  std::string sensor_id_;  // sensor_id
  std::vector<int> channels_;  // channel
  std::vector<float> tlapmean_;

  std::vector<double> biascoeffs_;
  bool newpc4pred_;  //  controls preconditioning due to sat-bias correction term
  bool adp_anglebc_;  //  logical to turn off or on the variational radiance angle bias correction
  bool emiss_bc_;  //  logical to turn off or on the emissivity predictor

  std::vector<std::string> predictors_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASRADIANCEGSI_H_
