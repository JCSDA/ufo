/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_
#define UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/marine/insitutemperature/ObsInsituTemperature.interface.h"
#include "ufo/operators/marine/insitutemperature/ObsInsituTemperatureParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// InsituTemperature observation operator class
class ObsInsituTemperature : public ObsOperatorBase,
                   private util::ObjectCounter<ObsInsituTemperature> {
 public:
  typedef ObsInsituTemperatureParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsInsituTemperature";}

  ObsInsituTemperature(const ioda::ObsSpace &, const ObsInsituTemperatureParameters &);
  ~ObsInsituTemperature() override;

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;


  const oops::Variables & requiredVars() const override {return varin_;}
  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  oops::ObsVariables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_
