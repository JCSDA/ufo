/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURETLAD_H_
#define UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURETLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/marine/insitutemperature/ObsInsituTemperatureParameters.h"
#include "ufo/operators/marine/insitutemperature/ObsInsituTemperatureTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// InsituTemperature for observation operator TL and AD class
class ObsInsituTemperatureTLAD : public LinearObsOperatorBase,
                                 private util::ObjectCounter<ObsInsituTemperatureTLAD> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsInsituTemperatureParameters Parameters_;
  static const std::string classname() {return "ufo::ObsInsituTemperatureTLAD";}

  ObsInsituTemperatureTLAD(const ioda::ObsSpace &, const ObsInsituTemperatureParameters &);
  ~ObsInsituTemperatureTLAD() override;

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}
  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
  oops::ObsVariables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURETLAD_H_
