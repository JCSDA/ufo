/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDENSITY_H_
#define UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDENSITY_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/OceanConversions/OceanConversions.interface.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling OceanDensity variable transform
class OceanDensityParameters : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(OceanDensityParameters, VariableTransformParametersBase)

 public:
  /// Input salinity
  oops::Parameter<std::string> SalinityVariable{"ocean salinity variable",
                                                "absoluteSalinity", this};
  oops::Parameter<std::string> SalinityGroup{"ocean salinity group",
                                             "ObsValue", this};
  /// Input temperature
  oops::Parameter<std::string> TemperatureVariable{"ocean temperature variable",
                                                   "waterTemperature", this};
  oops::Parameter<std::string> TemperatureGroup{"ocean temperature group",
                                                "ObsValue", this};
  /// Input pressure
  oops::Parameter<std::string> PressureVariable{"ocean pressure variable",
                                                "waterPressure", this};
  oops::Parameter<std::string> PressureGroup{"ocean pressure group",
                                             "ObsValue", this};
  /// Output density variable name
  oops::Parameter<std::string> DensityVariable{"ocean density variable",
                                               "waterDensity", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs density calculated from salinity, temperature, pressure.
///
/// Example
///
/// \code{.yaml}
/// obs filters:
/// - filter: Variable Transforms
///   Transform: OceanDensity
///   ocean pressure variable: waterPressure
///   ocean pressure group: DerivedObsValue
///   ocean temperature variable: waterTemperature
///   ocean temperature group: ObsValue
///   ocean salinity variable: absoluteSalinity
///   ocean salinity group: ObsValue
/// \endcode
///
/// will return density (kg/m^3) in a variable named (by default) "DerivedObsValue/waterDensity",
/// given salinity (g/kg), temperature (deg.C) and pressure (dbar), as given by TEOS-10
/// function gsw_rho_t_exact.
///

class OceanDensity : public TransformBase {
 public:
  typedef OceanDensityParameters Parameters_;
  OceanDensity(const Parameters_ &options,
               const ObsFilterData &data,
               const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
               const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string salinityvariable_;
  std::string salinitygroup_;
  std::string temperaturevariable_;
  std::string temperaturegroup_;
  std::string pressurevariable_;
  std::string pressuregroup_;
  std::string densityvariable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDENSITY_H_
