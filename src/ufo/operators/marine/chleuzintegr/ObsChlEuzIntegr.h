/*
 * (C) Copyright 2020-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_
#define UFO_OPERATORS_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"

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
// Chlorophyll Ocean Color Observation Operator
// Chlorophyll concentration (mg/m3) averaged over the euphotic layer (euz) that
// satellite sensors can see through
// euz is estimated based on surface chlorophyll concentration

class ObsChlEuzIntegrParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsChlEuzIntegrParameters, ObsOperatorParametersBase)
  // No additional options in YAML
};

class ObsChlEuzIntegr : public ObsOperatorBase,
                   private util::ObjectCounter<ObsChlEuzIntegr> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsChlEuzIntegr";}
  typedef ObsChlEuzIntegrParameters Parameters_;

  ObsChlEuzIntegr(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsChlEuzIntegr();

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
#endif  // UFO_OPERATORS_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_
