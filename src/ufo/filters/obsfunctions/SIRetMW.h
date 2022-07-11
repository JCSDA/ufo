/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SIRETMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_SIRETMW_H_

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
/// \brief Options applying to the retrieval of scattering index from 90.0 GHz and
//  150.0 GHz channels.
///
class SIRetMWParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SIRetMWParameters, Parameters)

 public:
  /// channel number corresponding to 90.0 GHz to which the retrieval
  /// of scattering index applies
  /// Example: MHS channel numbers used in scattering index retrieval
  ///          siret_channels: 1
  oops::RequiredParameter<int> ch90{"siret_ch90", this};

  /// channel number corresponding to 150.0 GHz to which the retrieval
  /// of scattering index applies
  /// Example: MHS channel numbers used in scattering index retrieval
  ///          siret_channels: 2
  oops::RequiredParameter<int> ch150{"siret_ch150", this};

  /// Names of the data group used to retrieve the scattering index
  /// Example: get retrieved SI from observation and simulated observation respectively
  ///          siret_types: [ObsValue, HofX]
  /// Example: get retrieved SI from observation or simulated observation only
  ///          siret_types: [ObsValue]
  ///          siret_types: [HofX]
  oops::RequiredParameter<std::vector<std::string>> varGroup{"siret_types", this};

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
/// \brief Retrieve scattering index from MHS 89.0 GHz and 157.0 GHz channels.
///
/// Reference: Geer et al. (2014)
/// Geer, A. J., Fabrizio, B., Bormann, N., & English, S. (2014). All-sky assimilation of
/// microwave humidity sounders. European Centre for Medium-Range Weather Forecasts.
///
class SIRetMW : public ObsFunctionBase<float> {
 public:
  explicit SIRetMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~SIRetMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  const std::vector<std::string> &siVariableGroups() const {
    return options_.varGroup.value();
  }
  inline static float getBadValue() {return bad_siret_value_;}
 private:
  ufo::Variables invars_;
  SIRetMWParameters options_;
  static constexpr float bad_siret_value_ = 1000.f;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SIRETMW_H_
