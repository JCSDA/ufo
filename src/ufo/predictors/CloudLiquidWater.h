/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_CLOUDLIQUIDWATER_H_
#define UFO_PREDICTORS_CLOUDLIQUIDWATER_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ObsBias;

// -----------------------------------------------------------------------------

/// Configuration parameters of the CloudLiquidWater predictor.
class CloudLiquidWaterParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(CloudLiquidWaterParameters, PredictorParametersBase)

 public:
    /// We must specify a sensor reference name such as SSMIS to know which channels to expect.
    oops::RequiredParameter<std::string> sensor{"sensor", this};
    /// In case we need to override the ObsValue group name with another optional group name.
    oops::Parameter<std::string> varGroup{"varGroup", "ObsValue", this};
    /// List below is solely for SSMIS data, but a different list of channel numbers could be
    /// added for a different sensor platform in the future.
    oops::OptionalParameter<int> ch19h{"ch19h", this};
    oops::OptionalParameter<int> ch19v{"ch19v", this};
    oops::OptionalParameter<int> ch22v{"ch22v", this};
    oops::OptionalParameter<int> ch37h{"ch37h", this};
    oops::OptionalParameter<int> ch37v{"ch37v", this};
    oops::OptionalParameter<int> ch91h{"ch91h", this};
    oops::OptionalParameter<int> ch91v{"ch91v", this};
    /// List below is solely for AMSU-A and ATMS data
    oops::OptionalParameter<int> ch238d{"clwdif_ch238", this};
    oops::OptionalParameter<int> ch314d{"clwdif_ch314", this};
    /// In case the inputs for 37v/h are clwret_ch37v and clwret_ch37h
    oops::OptionalParameter<int> clwret_ch37v{"clwret_ch37v", this};
    oops::OptionalParameter<int> clwret_ch37h{"clwret_ch37h", this};
    oops::OptionalParameter<int> order{"order", this};
    /// Path to an tlapse rate file for GMI only.
    oops::OptionalParameter<std::string> tlapse{"tlapse", this};
};

// -----------------------------------------------------------------------------

class CloudLiquidWater : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef CloudLiquidWaterParameters Parameters_;

  CloudLiquidWater(const Parameters_ &, const oops::ObsVariables &);
  ~CloudLiquidWater() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

  static void clwDerivative_amsua(const std::vector<float> &,
                           const std::vector<float> &,
                           const std::vector<float> &,
                           const std::vector<float> &,
                           const std::vector<float> &,
                           const std::vector<float> &,
                           std::vector<float> &);
  static void clw_bias_correction_gmi(const ObsBias &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<float> &,
                                      const std::vector<int> &,
                                      const std::vector<std::vector<std::vector<float>>> &,
                                      const std::vector<std::vector<float>> &,
                                      const std::vector<float> &,
                                      const std::vector<int> &,
                                      const  int &,
                                      std::vector<float> &,
                                      std::vector<float> &);

 private:
  CloudLiquidWaterParameters options_;
  std::vector<int> channels_;
  int order_;
  std::map<int, float> tlapmean_;  // <channel, tlaps>
  static constexpr float bad_clwret_value_ = 1000.f;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_CLOUDLIQUIDWATER_H_
