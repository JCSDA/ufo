/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOSTLAD_H_
#define UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOSTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/aerosols/MR/ObsAodGeos.h"
#include "ufo/operators/aerosols/MR/ObsAodGeosTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// AodGeos TL/AD observation operator class
class ObsAodGeosTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAodGeosTLAD> {
 public:
  typedef ObsAodGeosParameters Parameters_;
  static const std::string classname() {return "ufo::ObsAodGeosTLAD";}

  ObsAodGeosTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodGeosTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOSTLAD_H_
