/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_HEIGHTFROMPRESSURE_H_
#define UFO_VARIABLETRANSFORMS_CAL_HEIGHTFROMPRESSURE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the pressure to height conversion.
class Cal_HeightFromPressureParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_HeightFromPressureParameters, VariableTransformParametersBase);

 public:
  /// Height coordinate name.
  oops::RequiredParameter<std::string> heightCoord{"height coordinate", this};

  /// Height coordinate group.
  oops::RequiredParameter<std::string> heightGroup{"height group", this};

  /// Pressure coordinate name.
  oops::RequiredParameter<std::string> pressureCoord{"pressure coordinate", this};

  /// Pressure coordinate group.
  oops::RequiredParameter<std::string> pressureGroup{"pressure group", this};
};

/*!
* \brief Converts pressures to heights.
*/
class Cal_HeightFromPressure : public TransformBase {
 public:
  typedef Cal_HeightFromPressureParameters Parameters_;

  Cal_HeightFromPressure(const Parameters_ &options,
                         const ObsFilterData &data,
                         const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                         const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;

 private:
  /// Height coordinate name.
  std::string heightCoord_;

  /// Height group name.
  std::string heightGroup_;

  /// Pressure coordinate name.
  std::string pressureCoord_;

  /// Pressure group name.
  std::string pressureGroup_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_HEIGHTFROMPRESSURE_H_
