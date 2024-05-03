/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_H_
#define UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/aerosols/AODExt/ObsAodExtParameters.h"
#include "ufo/operators/aerosols/AODExt/ObsAodExtTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// AodExt TL/AD observation operator class
class ObsAodExtTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAodExtTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAodExtTLAD";}
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsAodExtParameters Parameters_;

  ObsAodExtTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodExtTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

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
#endif  // UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTTLAD_H_
