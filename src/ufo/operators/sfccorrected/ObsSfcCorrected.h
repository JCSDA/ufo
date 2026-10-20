/*
 * (C) Crown Copyright 2024, Met Office
 * (C) Copyright 2024 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTED_H_
#define UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTED_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/ObsOperatorBase.h"
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
/// SfcCorrected observation operator class
class ObsSfcCorrected : public ObsOperatorBase,
                   private util::ObjectCounter<ObsSfcCorrected> {
 public:
  typedef ObsSfcCorrectedParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsSfcCorrected";}

  ObsSfcCorrected(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsSfcCorrected();

  // Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

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

#endif  // UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTED_H_
