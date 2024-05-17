/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERP_H_
#define UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERP_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/atmsfcinterp/ObsAtmSfcInterp.interface.h"
#include "ufo/operators/atmsfcinterp/ObsAtmSfcInterpParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AtmSfcInterp observation operator class
class ObsAtmSfcInterp : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAtmSfcInterp> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsAtmSfcInterp";}
  typedef ObsAtmSfcInterpParameters Parameters_;

  ObsAtmSfcInterp(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAtmSfcInterp();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  oops::ObsVariables simulatedVars() const override {return operatorVars_;}

  int & toFortran() {return keyOperAtmSfcInterp_;}
  const int & toFortran() const {return keyOperAtmSfcInterp_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAtmSfcInterp_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  oops::ObsVariables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERP_H_
