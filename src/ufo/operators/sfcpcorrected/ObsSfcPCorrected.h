/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTED_H_
#define UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTED_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/sfcpcorrected/ObsSfcPCorrectedParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// SfcPCorrected observation operator class
class ObsSfcPCorrected : public ObsOperatorBase,
                   private util::ObjectCounter<ObsSfcPCorrected> {
 public:
  typedef ObsSfcPCorrectedParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsSfcPCorrected";}

  ObsSfcPCorrected(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsSfcPCorrected();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  oops::ObsVariables simulatedVars() const override {return operatorVars_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  oops::ObsVariables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTED_H_
