/*
 * (C) Copyright 2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AEROSOLS_AOP_OBSAODEXT_H_
#define UFO_AEROSOLS_AOP_OBSAODEXT_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/aerosols/AOP/ObsAodExt.interface.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AodExt observation operator class
class ObsAodExt : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAodExt> {
 public:
  static const std::string classname() {return "ufo::ObsAodExt";}

  ObsAodExt(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAodExt();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_AEROSOLS_AOP_OBSAODEXT_H_
