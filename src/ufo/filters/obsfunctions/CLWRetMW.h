/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to the retrieval of cloud liquid water from 23.8 GHz and
//  31.4 GHz channels.
///
class CLWRetMWParameters : public oops::Parameters {
 public:
  /// channel number corresponding to 23.8 GHz to which the retrieval
  /// of cloud liquid water applies
  /// Example: AMSU-A channel numbers used in cloud liquid water retrieval
  ///          clwret_channels: 1
  oops::RequiredParameter<int> ch238{"clwret_ch238", this};

  /// channel number corresponding to 31.4 GHz to which the retrieval
  /// of cloud liquid water applies
  /// Example: AMSU-A channel numbers used in cloud liquid water retrieval
  ///          clwret_channels: 2
  oops::RequiredParameter<int> ch314{"clwret_ch314", this};

  /// Names of the data group used to retrieve the cloud liquid water
  /// Example: get retrieved CLW from observation and simulated observation respectively
  ///          clwret_types: [ObsValue, HofX]
  /// Example: get retrieved CLW from observation or simulated observation only
  ///          clwret_types: [ObsValue]
  ///          clwret_types: [HofX]
  oops::RequiredParameter<std::vector<std::string>> varGroup{"clwret_types", this};

  /// Name of the data group to which the bias correction is applied (default is HofX)
  /// Example: add bias corretion to simulated observation
  ///          bias_application: HofX
  /// Example: add bias corretion to observation
  ///          bias_application: ObsValue
  oops::Parameter<std::string> addBias{"bias_application", "HofX", this};

  /// Name of the bias correction group used to replace the default group (default is ObsBias)
  /// Example: use observation bias correction values from GSI
  ///          test_bias: GsiObsBias
  oops::Parameter<std::string> testBias{"test_bias", "ObsBias", this};
};

///
/// \brief Retrieve cloud liquid water from 23.8 GHz and 31.4 GHz channels.
///
/// Reference: Grody et al. (2001)
/// Determination of precipitable water and cloud liquid water over oceans from
/// the NOAA 15 advanced microwave sounding unit
/// Journal of Geophysical Research (Vol. 106, No. D3, Pages 2943-2953)
///
class CLWRetMW : public ObsFunctionBase {
 public:
  explicit CLWRetMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~CLWRetMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  const std::vector<std::string> &clwVariableGroups() const {
    return options_.varGroup.value();
  }
  static void cloudLiquidWater(const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               std::vector<float> &,
                               const std::size_t);
  inline static float getBadValue() {return bad_clwret_value_;}
 private:
  ufo::Variables invars_;
  CLWRetMWParameters options_;
  static constexpr float bad_clwret_value_ = 1000.f;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_H_
