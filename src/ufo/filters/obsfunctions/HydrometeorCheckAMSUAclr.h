/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKAMSUACLR_H_
#define UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKAMSUACLR_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

///
/// \brief Options controling the cloud and precipitation checks for WM sensors
///
class HydrometeorCheckAMSUAclrParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HydrometeorCheckAMSUAclrParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the checks apply
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"test_hofx", "HofX", this};

  /// Name of the group for bias correction used to replace the default group (default is
  /// ObsBiasData)
  /// Example: use observation bias correction values from GSI
  ///          test_bias: GsiObsBias
  oops::Parameter<std::string> testBias{"test_bias", "ObsBiasData", this};

  /// Name of the group for bias correction terms used to replace the default group
  /// (default is ObsBiasTerm)
  /// Example: use observation bias correction terms from GSI
  ///          test_biasterm: GsiObsBiasTerm
  oops::Parameter<std::string> testBiasTerm{"test_biasterm", "ObsBiasTerm", this};

  /// Name of the group for bias correction predictors used to replace the default group
  /// (default is cloud_liquid_waterPredictor)
  oops::OptionalParameter<std::string> testBiasPred{"test_biaspredictor", this};
};

///
/// \brief Cloud and precipitation checks for AMSUA
/// Checks for all observations:
///   (1) Sanity check on CLW derivative values
/// Checks observation over all surface types include:
///   (1) Scattering check based on 54.4GHz channel
///   (2) Thick cloud check based on 52.8GHz channel
///   (3) Sensitivity to surface emissivity
/// Output of this function:
///    0 = channel is not affected by thick clouds and precipitation
///    1 = channel is affected by thick clouds and precipitataion
///
class HydrometeorCheckAMSUAclr : public ObsFunctionBase<float> {
 public:
  explicit HydrometeorCheckAMSUAclr(const eckit::LocalConfiguration &);
  ~HydrometeorCheckAMSUAclr();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  HydrometeorCheckAMSUAclrParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKAMSUACLR_H_
