/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_ADT_OBSADTTLAD_H_
#define UFO_OPERATORS_MARINE_ADT_OBSADTTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/marine/adt/ObsADT.h"
#include "ufo/operators/marine/adt/ObsADTTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// ADT for observation operator TL and AD class
class ObsADTTLAD : public LinearObsOperatorBase,
                   private util::ObjectCounter<ObsADTTLAD> {
 public:
  typedef ObsADTParameters Parameters_;
  static const std::string classname() {return "ufo::ObsADTTLAD";}

  ObsADTTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsADTTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_ADT_OBSADTTLAD_H_
