/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_
#define UFO_OPERATORS_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/marine/seaicethickness/ObsSeaIceThickness.interface.h"

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
class ObsSeaIceThicknessParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsSeaIceThicknessParameters, ObsOperatorParametersBase )
};
/// Sea ice thickness observation operator class
class ObsSeaIceThickness : public ObsOperatorBase,
                           private util::ObjectCounter<ObsSeaIceThickness> {
 public:
  typedef ObsSeaIceThicknessParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsSeaIceThickness";}

  ObsSeaIceThickness(const ioda::ObsSpace &, const ObsSeaIceThicknessParameters &);
  virtual ~ObsSeaIceThickness();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

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
#endif  // UFO_OPERATORS_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_
