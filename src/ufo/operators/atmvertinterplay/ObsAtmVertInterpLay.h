/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_H_
#define UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/atmvertinterplay/ObsAtmVertInterpLay.interface.h"
#include "ufo/operators/atmvertinterplay/ObsAtmVertInterpLayParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AtmVertInterpLay observation operator class
class ObsAtmVertInterpLay : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAtmVertInterpLay> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsAtmVertInterpLay";}

  typedef ObsAtmVertInterpLayParameters Parameters_;
  ObsAtmVertInterpLay(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAtmVertInterpLay();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAtmVertInterpLay_;}
  const int & toFortran() const {return keyOperAtmVertInterpLay_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAtmVertInterpLay_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAY_H_
