/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 - UK Met Office.
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
#include "ufo/operators/logarithm/ObsLogarithmParameters.h"

/// Forward declarations
namespace ioda {
class ObsSpace;
class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

/// \brief Logarithm observation operator.
///
/// This observation operator computes the H(x) vector by taking the logarithm
/// (with a chosen base) of the GeoVaLs variable at the lowest model level
/// (similar to the Identity operator). If the log base is not set, the natural
/// logarithm (base e) is used. If values are zero or negative, the output is
/// set to missing.
///
/// An example yaml configuration is:
///
///     obs operator:
///      name: Logarithm
///      base: 10
///
/// The associated parameters class has details of all available options.
class ObsLogarithm : public ObsOperatorBase,
                     private util::ObjectCounter<ObsLogarithm> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() { return "ufo::ObsLogarithm"; }
  typedef ObsLogarithmParameters Parameters_;

  ObsLogarithm(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsLogarithm() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
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

  /// Ref to the ObsSpace
  const ioda::ObsSpace &odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
