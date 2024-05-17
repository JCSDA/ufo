/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_VERTINTERP_OBSVERTINTERP_H_
#define UFO_OPERATORS_VERTINTERP_OBSVERTINTERP_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/vertinterp/ObsVertInterpParameters.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// VertInterp observation operator
class ObsVertInterp : public ObsOperatorBase,
                      private util::ObjectCounter<ObsVertInterp> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsVertInterp";}

  typedef ObsVertInterpParameters Parameters_;

  ObsVertInterp(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsVertInterp();

/// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  oops::ObsVariables simulatedVars() const override {return operatorVars_;}

  int & toFortran() {return keyOperVertInterp_;}
  const int & toFortran() const {return keyOperVertInterp_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperVertInterp_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  oops::ObsVariables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_VERTINTERP_OBSVERTINTERP_H_
