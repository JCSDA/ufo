/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_IDENTITY_OBSIDENTITYTLAD_H_
#define UFO_IDENTITY_OBSIDENTITYTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/identity/ObsIdentityTLAD.interface.h"
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
/// Identity TL/AD observation operator class
class ObsIdentityTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsIdentityTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsIdentityTLAD";}

  ObsIdentityTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsIdentityTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  oops::Variables simulatedVars() const override {return operatorVars_;}

  int & toFortran() {return keyOperObsIdentity_;}
  const int & toFortran() const {return keyOperObsIdentity_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperObsIdentity_;
  oops::Variables varin_;
  oops::Variables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_IDENTITY_OBSIDENTITYTLAD_H_
