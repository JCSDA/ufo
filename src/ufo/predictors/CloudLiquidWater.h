/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_CLOUDLIQUIDWATER_H_
#define UFO_PREDICTORS_CLOUDLIQUIDWATER_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

///
/// \brief Option to override varGroup default of ObsValue for all channels with ObsBias
///
class CloudLiquidWaterParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CloudLiquidWaterParameters, Parameters)

 public:
  /// We must specify a satellite reference name such as SSMIS to know which channels to expect.
  oops::RequiredParameter<std::string> satellite{"satellite", this};
  /// In case we need to override the ObsValue group name with another optional group name.
  oops::Parameter<std::string> varGroup{"varGroup", "ObsValue", this};
  /// List below is solely for SSMIS data, but a different list of channel numbers could be
  /// added for a different satellite platform in the future.
  oops::OptionalParameter<int> ch19h{"ch19h", this};
  oops::OptionalParameter<int> ch19v{"ch19v", this};
  oops::OptionalParameter<int> ch22v{"ch22v", this};
  oops::OptionalParameter<int> ch37h{"ch37h", this};
  oops::OptionalParameter<int> ch37v{"ch37v", this};
  oops::OptionalParameter<int> ch91h{"ch91h", this};
  oops::OptionalParameter<int> ch91v{"ch91v", this};
};

// -----------------------------------------------------------------------------

class CloudLiquidWater : public PredictorBase {
 public:
  CloudLiquidWater(const eckit::Configuration &, const oops::Variables &);
  ~CloudLiquidWater() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  CloudLiquidWaterParameters options_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_CLOUDLIQUIDWATER_H_
