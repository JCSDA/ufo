/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_PROFILEHORIZONTALDRIFT_H_
#define UFO_VARIABLETRANSFORMS_CAL_PROFILEHORIZONTALDRIFT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/*!
* \brief Profile horizontal drift calculation.
*
* This computes the horizontal drift positions (and times) given the
* horizontal wind speeds, heights, and an assumed rate of ascent.
*
* This function should only be applied to sondes whose horizontal position was not already
* measured (i.e. not to sondes reporting in BUFR format).
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transform
*   Transform: ["ProfileHorizontalDrift"]
* \endcode
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_ProfileHorizontalDrift : public TransformBase {
 public:
  Cal_ProfileHorizontalDrift(const GenericVariableTransformParameters &options,
                             const ObsFilterData &data,
                             const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_PROFILEHORIZONTALDRIFT_H_
