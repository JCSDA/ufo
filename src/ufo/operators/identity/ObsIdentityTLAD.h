/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_IDENTITY_OBSIDENTITYTLAD_H_
#define UFO_OPERATORS_IDENTITY_OBSIDENTITYTLAD_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/identity/ObsIdentityParameters.h"

// Forward declarations
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
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() { return "ufo::ObsIdentityTLAD"; }
  typedef ObsIdentityParameters Parameters_;

  ObsIdentityTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsIdentityTLAD();

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables &requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  /// Level index 0 is closest to surface.
  bool levelIndexZeroAtSurface_;
};


}  // namespace ufo
#endif  // UFO_OPERATORS_IDENTITY_OBSIDENTITYTLAD_H_
