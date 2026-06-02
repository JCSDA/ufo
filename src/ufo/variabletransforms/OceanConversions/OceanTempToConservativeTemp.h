/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANTEMPTOCONSERVATIVETEMP_H_
#define UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANTEMPTOCONSERVATIVETEMP_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/OceanConversions/OceanConversions.interface.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling OceanTempToConservativeTemp variable transform
class OceanTempToConservativeTempParameters : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(OceanTempToConservativeTempParameters, VariableTransformParametersBase)

 public:
  /// Input salinity
  oops::Parameter<std::string> SalinityVariable{"ocean salinity variable",
                                                "absoluteSalinity", this};
  oops::Parameter<std::string> SalinityGroup{"ocean salinity group",
                                             "ObsValue", this};
  /// Input in-situ temperature
  oops::Parameter<std::string> TemperatureVariable{"ocean temperature variable",
                                                   "waterTemperature", this};
  oops::Parameter<std::string> TemperatureGroup{"ocean temperature group",
                                                "ObsValue", this};
  /// Input pressure
  oops::Parameter<std::string> PressureVariable{"ocean pressure variable",
                                                "waterPressure", this};
  oops::Parameter<std::string> PressureGroup{"ocean pressure group",
                                             "ObsValue", this};
  /// Output water conservative temperature variable name
  oops::Parameter<std::string> ConservativeTempVariable{"ocean conservative temperature name",
                                                        "waterConservativeTemperature", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs conservative temperature converted from in situ temperature,
/// using TEOS-10 function gsw_ct_from_t.
///
/// Example
///
/// \code{.yaml}
/// obs filters:
/// - filter: Variable Transforms
///   Transform: OceanTempToConservativeTemp
///   ocean pressure variable: waterPressure
///   ocean pressure group: DerivedObsValue
///   ocean temperature variable: waterTemperature
///   ocean temperature group: ObsValue
///   ocean salinity variable: salinity
///   ocean salinity group: ObsValue
/// \endcode
///
/// will return conservative temperature (deg.C) in a variable named (by default)
/// "DerivedObsValue/waterConservativeTemperature",
/// given absolute salinity (g/kg), in situ temperature (deg.C) and pressure (dbar).
///

class OceanTempToConservativeTemp : public TransformBase {
 public:
  typedef OceanTempToConservativeTempParameters Parameters_;
  OceanTempToConservativeTemp(const Parameters_ &options,
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
  std::string conservativetempvariable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANTEMPTOCONSERVATIVETEMP_H_
