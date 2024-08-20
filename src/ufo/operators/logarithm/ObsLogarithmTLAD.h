/*
 * (C) Copyright 2023- UCAR
 # (C) Crown Copyright 2024 - UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/logarithm/ObsLogarithmParameters.h"

// Forward declarations
namespace ioda {
class ObsSpace;
class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// TL/AD code for the Logarithm observation operator.
class ObsLogarithmTLAD : public LinearObsOperatorBase,
                         private util::ObjectCounter<ObsLogarithmTLAD> {
 public:
  static const std::string classname() { return "ufo::ObsLogarithmTLAD"; }
  typedef ObsLogarithmParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  ObsLogarithmTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsLogarithmTLAD();

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &,
                     const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &,
                     const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &,
                     const QCFlags_t &) const override;

  const oops::Variables &requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

 private:
  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  /// Log base (0 for natural logarithm)
  float logBase_ = 0.0;

  /// Values of derivative at the lowest model level - length of the vector is
  /// the number of required variables in the operator multiplied by the number
  /// of locations.
  std::vector<double> dHdx_;

  /// Ref to the ObsSpace
  const ioda::ObsSpace &odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
