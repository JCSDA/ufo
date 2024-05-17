/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_ADT_OBSADT_H_
#define UFO_OPERATORS_MARINE_ADT_OBSADT_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/marine/adt/ObsADTParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// ADT observation operator class
class ObsADT : public ObsOperatorBase,
                   private util::ObjectCounter<ObsADT> {
 public:
  typedef ObsADTParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsADT";}

  ObsADT(const ioda::ObsSpace &, const ObsADTParameters &);
  virtual ~ObsADT();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override {return requiredVars_;}
  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

  oops::Variables requiredVars_;
  const ioda::ObsSpace& odb_;
  oops::ObsVariables operatorVars_;
  int operatorVarIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_ADT_OBSADT_H_
