/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOS_H_
#define UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOS_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/aerosols/MR/ObsAodGeos.interface.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class ObsAodGeosParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAodGeosParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<std::vector<std::string>> tracerGeovals{
        "tracer_geovals", this};
  oops::RequiredParameter<std::vector<double>> wavelengths{
        "wavelengths", this};
  oops::RequiredParameter<std::string> rcfile{"RCFile", this};
};

// -----------------------------------------------------------------------------
/// AodGeos observation operator class
class ObsAodGeos : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAodGeos> {
 public:
  typedef ObsAodGeosParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAodGeos";}

  ObsAodGeos(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodGeos();

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
#endif  // UFO_OPERATORS_AEROSOLS_MR_OBSAODGEOS_H_
