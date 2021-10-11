/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RTTOV_OBSRADIANCERTTOV_H_
#define UFO_RTTOV_OBSRADIANCERTTOV_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/rttov/ObsRadianceRTTOV.interface.h"
#include "ufo/rttov/ObsRadianceRTTOVParameters.h"

namespace eckit {
  class Configuration;
  class LocalConfiguration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceRTTOV observation for UFO.
class ObsRadianceRTTOV : public ObsOperatorBase,
                    private util::ObjectCounter<ObsRadianceRTTOV> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceRTTOV";}

  ObsRadianceRTTOV(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceRTTOV();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperRadianceRTTOV_;}
  const int & toFortran() const {return keyOperRadianceRTTOV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperRadianceRTTOV_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsRadianceRTTOVParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RTTOV_OBSRADIANCERTTOV_H_
