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
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/product/ObsProductParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Product observation operator.
///
/// This observation operator computes the H(x) vector by multiplying the model values at the lowest
/// model level with another customizable GeoVaLs variable. Note that this operator acts similarly
/// to the Identity operator and only acts on the lowest model level. The scaling GeoVaLs must be
/// two dimensional to match the rank after choosing only the lowest model level.
///
/// An example yaml configuration is:
///
///     obs operator:
///      name: Product
///
/// The associated parameters class has details of all available options.
class ObsProduct : public ObsOperatorBase,
  private util::ObjectCounter<ObsProduct> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsProduct";}
  typedef ObsProductParameters Parameters_;

  ObsProduct(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsProduct() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

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

  /// Geoval name to act on
  std::string geovalName_ = "";

  /// Exponent of geovals
  float scalingVariableExponent_ = 0;

  /// Group for variable to scale hofx by
  std::string variableGroupToScaleHofxBy_;

  /// Variable to scale hofx by
  std::string variableNameToScaleHofxBy_;

  /// Ref to the ObsSpace
  const ioda::ObsSpace& odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
