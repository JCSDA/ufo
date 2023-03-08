/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANPRACTICALSALINITYTOABSOLUTESALINITY_H_
#define UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANPRACTICALSALINITYTOABSOLUTESALINITY_H_

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

/// \brief Options controlling OceanPracticalSalinityToAbsoluteSalinity variable transform
class OceanPracticalSalinityToAbsoluteSalinityParameters : public VariableTransformParametersBase
{
  OOPS_CONCRETE_PARAMETERS(OceanPracticalSalinityToAbsoluteSalinityParameters,
                           VariableTransformParametersBase)

 public:
  /// Input salinity
  oops::Parameter<std::string> PracticalSalinityVariable{"ocean practical salinity variable",
                                                         "salinity", this};
  oops::Parameter<std::string> PracticalSalinityGroup{"ocean practical salinity group",
                                                      "ObsValue", this};
  /// Input pressure
  oops::Parameter<std::string> PressureVariable{"ocean pressure variable",
                                                "waterPressure", this};
  oops::Parameter<std::string> PressureGroup{"ocean pressure group",
                                             "ObsValue", this};
  /// Output absolute salinity variable name
  oops::Parameter<std::string> AbsoluteSalinityVariable{"ocean absolute salinity name",
                                                        "absoluteSalinity", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs absolute salinity converted from practical salinity,
/// using TEOS-10 function gsw_SA_from_PS.
///
/// Example
///
/// \code{.yaml}
/// obs filters:
/// - filter: Variable Transforms
///   Transform: OceanPracticalSalinityToAbsoluteSalinity
///   ocean pressure variable: waterPressure
///   ocean pressure group: DerivedObsValue
///   ocean practical salinity variable: salinity
///   ocean practical salinity group: ObsValue
///   ocean absolute salinity variable: absoluteSalinity
/// \endcode
///
/// will return absolute salinity (g/kg) in a variable named (by default)
/// "DerivedObsValue/absoluteSalinity",
/// given practical salinity and pressure (dbar).
///

class OceanPracticalSalinityToAbsoluteSalinity : public TransformBase {
 public:
  typedef OceanPracticalSalinityToAbsoluteSalinityParameters Parameters_;
  OceanPracticalSalinityToAbsoluteSalinity(const Parameters_ &options,
                                     const ObsFilterData &data,
                                     const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                     const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string practicalsalinityvariable_;
  std::string practicalsalinitygroup_;
  std::string pressurevariable_;
  std::string pressuregroup_;
  std::string absolutesalinityvariable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANPRACTICALSALINITYTOABSOLUTESALINITY_H_
