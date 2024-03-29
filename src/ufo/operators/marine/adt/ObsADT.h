/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_MARINE_ADT_OBSADT_H_
#define UFO_OPERATORS_MARINE_ADT_OBSADT_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/operators/marine/adt/ObsADT.interface.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
class ObsADTParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsADTParameters, ObsOperatorParametersBase )
  // NO extra parameters needed
};

/// ADT observation operator class
class ObsADT : public ObsOperatorBase,
                   private util::ObjectCounter<ObsADT> {
 public:
  typedef ObsADTParameters Parameters_;
  static const std::string classname() {return "ufo::ObsADT";}

  ObsADT(const ioda::ObsSpace &, const ObsADTParameters &);
  virtual ~ObsADT();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_ADT_OBSADT_H_
