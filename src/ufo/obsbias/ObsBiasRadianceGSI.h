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

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
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
  class ObsBiasIncrement;
  class ObsDiagnostics;

/// Class to handle observation bias model from GSI Radiance.

class ObsBiasRadianceGSI : public ObsBiasBase,
                           private util::ObjectCounter<ObsBiasRadianceGSI> {
 public:
  static const std::string classname() {return "ufo::ObsBiasRadianceGSI";}

/// Constructor
  ObsBiasRadianceGSI(const ioda::ObsSpace &, const eckit::Configuration &);

/// Destructor
  virtual ~ObsBiasRadianceGSI() {}

/// I/O and diagnostics
  void read(const eckit::Configuration &) override;
  void write(const eckit::Configuration &) const override;
  double norm() const override;

/// Add increments
  ObsBiasRadianceGSI & operator+=(const ObsBiasIncrement &) override;
  ObsBiasRadianceGSI & operator=(const ObsBias &) override;

/// Obs bias operator
  void computeObsBias(ioda::ObsVector &,
                      ioda::ObsDataVector<float> &) const override;

/// Obs bias predictor
  void computeObsBiasPredictors(const GeoVaLs &, const ObsDiagnostics &,
                                ioda::ObsDataVector<float> &) const override;

/// Other
  const oops::Variables & requiredGeoVaLs() const override {return *geovars_;}
  const oops::Variables & requiredHdiagnostics() const override {return *hdiags_;}
  const oops::Variables & predNames() const override {return *predNames_;}

/// Bias parameters interface
  std::size_t size() const override {return biascoeffs_.size();}
  double & operator[](const unsigned int ii) override {return biascoeffs_[ii];}

  const ioda::ObsSpace & obspace() const override {return odb_;}
 private:
  void print(std::ostream &) const override;

  const ioda::ObsSpace & odb_;
  std::unique_ptr<const oops::Variables> geovars_;
  std::unique_ptr<const oops::Variables> hdiags_;
  std::unique_ptr<const oops::Variables> predNames_;
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
