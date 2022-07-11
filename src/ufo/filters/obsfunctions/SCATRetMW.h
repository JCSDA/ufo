/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SCATRETMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_SCATRETMW_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options applying to the retrieval of scattering index from 23.8 GHz, 31.4 GHz,
//  and 89 GHz channels.
///
class SCATRetMWParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SCATRetMWParameters, Parameters)

 public:
  /// channel number corresponding to 23.8 GHz to which the retrieval
  /// of scattering index applies
  /// Example: AMSU-A channel numbers used in scattering index retrieval
  ///          scatret_channels: 1
  oops::RequiredParameter<int> ch238{"scatret_ch238", this};

  /// channel number corresponding to 31.4 GHz to which the retrieval
  /// of scattering index applies
  /// Example: AMSU-A channel numbers used in scattering index retrieval
  ///          scatret_channels: 2
  oops::RequiredParameter<int> ch314{"scatret_ch314", this};

  /// channel number corresponding to 89 GHz to which the retrieval
  /// of scattering index applies
  /// Example: AMSU-A channel numbers used in scattering index retrieval
  ///          scatret_channels: 15
  oops::RequiredParameter<int> ch890{"scatret_ch890", this};

  /// Names of the data group used to retrieve the scattering index
  /// Example: get retrieved scattering index from observation and simulated values respectively
  ///          scatret_types: [ObsValue, HofX]
  /// Example: get retrieved scattering index from observation or simulated values only
  ///          scatret_types: [ObsValue]
  ///          scatret_types: [HofX]
  oops::RequiredParameter<std::vector<std::string>> varGroup{"scatret_types", this};

  /// Name of the data group to which the bias correction is applied (default is HofX)
  /// Example: add bias corretion to simulated observation
  ///          bias_application: HofX
  /// Example: add bias corretion to observation
  ///          bias_application: ObsValue
  oops::Parameter<std::string> addBias{"bias_application", "HofX", this};

  /// Name of the bias correction group used to replace the default group (default is ObsBiasData)
  /// Example: use observation bias correction values from GSI
  ///          test_bias: GsiObsBias
  oops::Parameter<std::string> testBias{"test_bias", "ObsBiasData", this};
};

///
/// \brief Retrieve scattering index from 23.8 GHz, 31.4 GHz and 89 GHz channels.
///
/// Reference: Grody et al. (2000)
/// Application of AMSU for obtaining hydrological parameters
/// Microw. Radiomet. Remote Sens. Eatch's Surf. Atmosphere, pp. 339-351
///
class SCATRetMW : public ObsFunctionBase<float> {
 public:
  explicit SCATRetMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~SCATRetMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  const std::vector<std::string> &scatVariableGroups() const {
    return options_.varGroup.value();
  }
  inline static float getBadValue() {return bad_scatret_value_;}
 private:
  ufo::Variables invars_;
  SCATRetMWParameters options_;
  static constexpr float bad_scatret_value_ = 1000.f;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SCATRETMW_H_
