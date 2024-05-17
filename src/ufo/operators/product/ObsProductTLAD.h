/*
 * (C) Copyright 2023- UCAR
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
#include "ufo/operators/product/ObsProductParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// TL/AD code for the Product observation operator.
class ObsProductTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsProductTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsProductTLAD";}
  typedef ObsProductParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  ObsProductTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsProductTLAD();

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override {return requiredVars_;}

  oops::ObsVariables simulatedVars() const override {return operatorVars_;}

 private:
  void print(std::ostream &) const override;

 private:
  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  /// Geoval name to act on
  std::string geovalName_ = "";

  /// Scaling variable
  std::vector<double> scalingVariable_;

  /// Group for scaling variable
  std::string variableGroupToScaleHofxBy_;

  /// Name for scaling variable
  std::string variableNameToScaleHofxBy_;

  /// Exponent of scaling variable
  float scalingVariableExponent_ = 0;

  /// Ref to the ObsSpace
  const ioda::ObsSpace& odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
