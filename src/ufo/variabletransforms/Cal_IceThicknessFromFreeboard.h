/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_ICETHICKNESSFROMFREEBOARD_H_
#define UFO_VARIABLETRANSFORMS_CAL_ICETHICKNESSFROMFREEBOARD_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling Cal_IceThicknessFromFreeboard variable transform
class Cal_IceThicknessFromFreeboardParameters
    : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_IceThicknessFromFreeboardParameters,
                           VariableTransformParametersBase)

 public:
  /// Input ice Freeboard
  oops::Parameter<Variable> IceFreeboardVariable{
      "ice freeboard variable", Variable("ObsValue/seaIceFreeboard"), this};
  /// Input sea water density
  oops::Parameter<Variable> WaterDensityVariable{
      "sea water density variable", Variable("ObsValue/seaWaterDensity"), this};
  /// Input snow depth
  oops::Parameter<Variable> SnowDepthVariable{
      "snow depth variable", Variable("ObsValue/totalSnowDepth"), this};
  /// Input snow density
  oops::Parameter<Variable> SnowDensityVariable{
      "snow density variable", Variable("ObsValue/snowDensity"), this};
  /// Input ice density
  oops::Parameter<Variable> IceDensityVariable{
      "ice density variable", Variable("ObsValue/iceDensity"), this};
  /// Output ice thickness variable name
  oops::Parameter<Variable> IceThicknessVariable{
      "ice thickness variable", Variable("iceThickness"), this};
  /// Output ice thickness error standard deviations
  oops::Parameter<bool> calculateErrors{"calculate error standard deviations",
                                        true, this};
  /// Input sea ice freeboard error standard deviation
  oops::Parameter<Variable> IceFreeboardESDVariable{
      "ice freeboard error standard deviation variable",
      Variable("ObservedErrorStandardDeviation/seaIceFreeboard"), this};
  /// Input ice density error standard deviation
  oops::Parameter<Variable> IceDensityESDVariable{
      "ice density error standard deviation variable",
      Variable("ObservedErrorStandardDeviation/iceDensity"), this};
  /// Input snow density error standard deviation
  oops::Parameter<Variable> SnowDensityESDVariable{
      "snow density error standard deviation variable",
      Variable("ObservedErrorStandardDeviation/snowDensity"), this};
  /// Input snow depth error standard deviation
  oops::Parameter<Variable> SnowDepthESDVariable{
      "snow depth error standard deviation variable",
      Variable("ObservedErrorStandardDeviation/totalSnowDepth"), this};
  /// Output ice thickness systematic error standard deviation
  oops::Parameter<Variable> IceThicknessSystematicESDVariable{
      "ice thickness systematic error standard deviation variable",
      Variable("SystematicErrorStandardDeviation/iceThickness"), this};
  /// Output ice thickness random error standard deviation
  oops::Parameter<Variable> IceThicknessRandomESDVariable{
      "ice thickness random error standard deviation variable",
      Variable("RandomErrorStandardDeviation/iceThickness"), this};
};

// -----------------------------------------------------------------------------

/// \brief Calculate ice thickness from ice freeboard, ice depth, and the
/// densities of the snow, ice, and surface water.
/// \details The ice thickness is calculated using the formula from
/// https://tc.copernicus.org/articles/8/1607/2014/tc-8-1607-2014.html:
///  ice thickness = (freeboard * water density + snow depth * snow density) /
///                  (water density - ice density)
///  ice thickness random esd = sqrt(
///    ( freeboard esd  * water density / (water density - ice density) )**2 +
///    ( ice density esd * (freeboard * water density + snow depth * snow
///    density) /
///      (water density - ice density)**2
///    )**2
///  )
///  ice thickness systematic esd = sqrt(
///    ( snow depth esd * snow density / (water density - ice density) )**2 +
///    ( snow density esd * snow depth / (water density - ice density) )**2
///  )
///
/// This will return ice thickness (m) in a variable named (by default)
/// "DerivedObsValue/iceThickness", given the following input variables:
/// - ice freeboard (m)
/// - water density (kg/m^3)
/// - snow depth (m)
/// - snow density (kg/m^3)
/// - ice density (kg/m^3)
/// The systematic and random errors of ice thickness are also calculated via
/// gaussian error propagation.
///
/// Example yaml:
///
/// \code{.yaml}
/// obs filters:
/// - filter: Variable Transforms
///   Transform: Calculate iceThickness from seaIceFreeboard
///   ice freeboard variable: ObsValue/seaIceFreeboard
///   sea water density variable: ObsValue/seaWaterDensity
///   snow depth variable: ObsValue/totalSnowDepth
///   snow density variable: ObsValue/snowDensity
///   ice density variable: ObsValue/iceDensity
///   ice freeboard error standard deviation variable:
///       ObservedErrorStandardDeviation/seaIceFreeboard
///   ice density error standard deviation variable:
///       ObservedErrorStandardDeviation/iceDensity
///   ice thickness variable: iceThickness
///   ice thickness systematic error standard deviation variable:
///       SystematicStandardDeviation/iceThickness
///   ice thickness random error standard deviation variable:
///       RandomStandardDeviation/iceThickness
/// \endcode
///

class Cal_IceThicknessFromFreeboard : public TransformBase {
 public:
  typedef Cal_IceThicknessFromFreeboardParameters Parameters_;
  Cal_IceThicknessFromFreeboard(
      const Parameters_ &options, const ObsFilterData &data,
      const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
      const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  const std::string classname_ = "Cal_IceThicknessFromFreeboard";
  Parameters_ options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_ICETHICKNESSFROMFREEBOARD_H_
