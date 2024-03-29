/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_MARINEVERTINTERP_OBSMARINEVERTINTERPTLAD_H_
#define UFO_OPERATORS_MARINE_MARINEVERTINTERP_OBSMARINEVERTINTERPTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/marine/marinevertinterp/ObsMarineVertInterpParameters.h"
#include "ufo/operators/marine/marinevertinterp/ObsMarineVertInterpTLAD.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// Marinevertinterp for observation operator TL and AD class
class ObsMarineVertInterpTLAD : public LinearObsOperatorBase,
                                 private util::ObjectCounter<ObsMarineVertInterpTLAD> {
 public:
  typedef ObsMarineVertInterpParameters Parameters_;
  static const std::string classname() {return "ufo::ObsMarineVertInterpTLAD";}

  ObsMarineVertInterpTLAD(const ioda::ObsSpace &, const ObsMarineVertInterpParameters &);
  ~ObsMarineVertInterpTLAD() override;

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}
  oops::Variables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
  oops::Variables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_MARINEVERTINTERP_OBSMARINEVERTINTERPTLAD_H_
