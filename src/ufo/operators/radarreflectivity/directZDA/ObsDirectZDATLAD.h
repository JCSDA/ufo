/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDATLAD_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDATLAD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/radarreflectivity/directZDA/ObsDirectZDA.h"
#include "ufo/operators/radarreflectivity/directZDA/ObsDirectZDATLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// DirectZDA TL/AD observation operator class
class ObsDirectZDATLAD : public LinearObsOperatorBase,
                         private util::ObjectCounter<ObsDirectZDATLAD> {
 public:
  typedef ObsDirectZDAParameters Parameters_;

  static const std::string classname() {return "ufo::ObsDirectZDATLAD";}

  ObsDirectZDATLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsDirectZDATLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperDirectZDA_;}
  const int & toFortran() const {return keyOperDirectZDA_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperDirectZDA_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDATLAD_H_
