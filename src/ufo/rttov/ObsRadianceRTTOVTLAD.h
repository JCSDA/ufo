/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_RTTOV_OBSRADIANCERTTOVTLAD_H_
#define UFO_RTTOV_OBSRADIANCERTTOVTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/rttov/ObsRadianceRTTOVTLAD.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceRTTOV TL/AD observation operator class
class ObsRadianceRTTOVTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsRadianceRTTOVTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceRTTOVTLAD";}

  ObsRadianceRTTOVTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceRTTOVTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperRadianceRTTOV_;}
  const int & toFortran() const {return keyOperRadianceRTTOV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperRadianceRTTOV_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RTTOV_OBSRADIANCERTTOVTLAD_H_
