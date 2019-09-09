/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITYTLAD_H_
#define UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITYTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/radarradialvelocity/ObsRadarRadialVelocityTLAD.interface.h"

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
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// RadarRadialVelocity TL/AD observation operator class
class ObsRadarRadialVelocityTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsRadarRadialVelocityTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsRadarRadialVelocityTLAD";}

  ObsRadarRadialVelocityTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadarRadialVelocityTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const;

  // Other
  const oops::Variables & variables() const {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITYTLAD_H_
