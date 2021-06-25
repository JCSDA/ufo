/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_IDENTITY_OBSIDENTITYTLAD_H_
#define UFO_IDENTITY_OBSIDENTITYTLAD_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

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
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// TL/AD code for the Identity observation operator.
class ObsIdentityTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsIdentityTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsIdentityTLAD";}

  ObsIdentityTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsIdentityTLAD();

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  const oops::Variables & requiredVars() const override {return requiredVars_;}

  oops::Variables simulatedVars() const override {return operatorVars_;}

 private:
  void print(std::ostream &) const override;

 private:
  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::Variables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_IDENTITY_OBSIDENTITYTLAD_H_
