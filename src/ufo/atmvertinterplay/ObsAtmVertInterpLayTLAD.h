/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYTLAD_H_
#define UFO_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmvertinterplay/ObsAtmVertInterpLayTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

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
/// AtmVertInterpLay observation operator
class ObsAtmVertInterpLayTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsAtmVertInterpLayTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAtmVertInterpLayTLAD";}

  ObsAtmVertInterpLayTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAtmVertInterpLayTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAtmVertInterpLay_;}
  const int & toFortran() const {return keyOperAtmVertInterpLay_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAtmVertInterpLay_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYTLAD_H_
