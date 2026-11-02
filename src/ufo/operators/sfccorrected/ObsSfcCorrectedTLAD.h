/*
 * (C) Crown Copyright 2024, Met Office
 * (C) Copyright 2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDTLAD_H_
#define UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/sfccorrected/ObsSfcCorrectedParameters.h"
#include "ufo/operators/sfccorrected/SurfaceOperatorBase.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// SfcCorrected observation Linear operator class
class ObsSfcCorrectedTLAD : public LinearObsOperatorBase,
                   private util::ObjectCounter<ObsSfcCorrectedTLAD> {
 public:
  using Parameters_ = ObsSfcCorrectedParameters;
  using QCFlags_t = ioda::ObsDataVector<int>;

  static std::string classname() {return "ufo::ObsSfcCorrectedTLAD";}

  ObsSfcCorrectedTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsSfcCorrectedTLAD();

  // Linear Obs Operator
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return requiredVars_;}
  oops::ObsVariables simulatedVars() const override {return operatorVars_;}

 private:
  void print(std::ostream &) const override;

  // Required variables from the GeoVaLs
  oops::Variables requiredVars_;

  // Link to odb
  const ioda::ObsSpace & odb_;

  // Parameters
  Parameters_ params_;

  // Operator variables
  oops::ObsVariables operatorVars_;

  // Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  // Vector of operators
  std::vector<std::unique_ptr<SurfaceOperatorBase>> operators_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDTLAD_H_
