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

/*!
* \brief Converts pressures to heights.
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_HeightFromPressure : public TransformBase {
 public:
  Cal_HeightFromPressure(const GenericVariableTransformParameters &options,
                         const ObsFilterData &data,
                         const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_HEIGHTFROMPRESSURE_H_
