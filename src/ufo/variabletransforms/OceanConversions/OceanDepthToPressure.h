/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDEPTHTOPRESSURE_H_
#define UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDEPTHTOPRESSURE_H_

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

/// \brief Options controlling OceanDepthToPressure variable transform
class OceanDepthToPressureParameters : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(OceanDepthToPressureParameters, VariableTransformParametersBase)

 public:
  /// Input depth variable of the ocean depth to pressure conversion
  oops::Parameter<std::string> DepthVariable{"ocean depth variable",
                                             "depthBelowWaterSurface", this};
  oops::Parameter<std::string> DepthGroup{"ocean depth group",
                                          "MetaData", this};
  /// Output pressure variable name
  oops::Parameter<std::string> PressureVariable{"ocean pressure name",
                                               "waterPressure", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs pressure converted from depth, using TEOS-10 function gsw_p_from_z.
///
/// Example
///
/// \code{.yaml}
/// obs filters:
/// - filter: Variable Transforms
///   Transform: OceanDepthToPressure
///   ocean depth variable: depthBelowWaterSurface
///   ocean depth group: DerivedObsValue
/// \endcode
///
/// will return pressure (dbar) in a variable named (by default) "DerivedObsValue/waterPressure",
/// as a function of depth (m) and latitude (deg).
/// Computed by TEOS-10 function gsw_p_from_z.
///

class OceanDepthToPressure : public TransformBase {
 public:
  typedef OceanDepthToPressureParameters Parameters_;
  OceanDepthToPressure(const Parameters_ &options,
                       const ObsFilterData &data,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);

    // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string depthvariable_;
  std::string depthgroup_;
  std::string pressurevariable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_OCEANCONVERSIONS_OCEANDEPTHTOPRESSURE_H_
