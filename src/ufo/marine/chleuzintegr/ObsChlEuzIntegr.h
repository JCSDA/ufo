/*
 * (C) Copyright 2020-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_
#define UFO_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

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
// Chlorophyll Ocean Color Observation Operator
// Chlorophyll concentration (mg/m3) averaged over the euphotic layer (euz) that
// satellite sensors can see through
// euz is estimated based on surface chlorophyll concentration

class ObsChlEuzIntegr : public ObsOperatorBase,
                   private util::ObjectCounter<ObsChlEuzIntegr> {
 public:
  static const std::string classname() {return "ufo::ObsChlEuzIntegr";}

  ObsChlEuzIntegr(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsChlEuzIntegr();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_CHLEUZINTEGR_OBSCHLEUZINTEGR_H_
