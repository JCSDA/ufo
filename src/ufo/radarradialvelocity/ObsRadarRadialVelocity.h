/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_
#define UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/radarradialvelocity/ObsRadarRadialVelocity.interface.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadarRadialVelocity observation operator class
class ObsRadarRadialVelocity : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadarRadialVelocity> {
 public:
  static const std::string classname() {return "ufo::ObsRadarRadialVelocity";}

  ObsRadarRadialVelocity(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadarRadialVelocity();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_
