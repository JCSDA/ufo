/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_SEAICEFRACTION_OBSSEAICEFRACTION_H_
#define UFO_OPERATORS_MARINE_SEAICEFRACTION_OBSSEAICEFRACTION_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorParametersBase.h"

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
class ObsSeaIceFractionParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsSeaIceFractionParameters, ObsOperatorParametersBase )
  // NO extra parameters needed
};

/// Total ice concentration observation operator class
class ObsSeaIceFraction : public ObsOperatorBase,
                          private util::ObjectCounter<ObsSeaIceFraction> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  typedef ObsSeaIceFractionParameters Parameters_;

  static const std::string classname() {return "ufo::ObsSeaIceFraction";}

  ObsSeaIceFraction(const ioda::ObsSpace &, const ObsSeaIceFractionParameters &);
  virtual ~ObsSeaIceFraction();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_SEAICEFRACTION_OBSSEAICEFRACTION_H_
