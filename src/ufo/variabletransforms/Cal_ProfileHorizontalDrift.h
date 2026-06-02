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

/// Configuration parameters for the profile horizontal drift calculation.
class Cal_ProfileHorizontalDriftParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_ProfileHorizontalDriftParameters, VariableTransformParametersBase);

 public:
  /// Height coordinate name.
  oops::RequiredParameter<std::string> HeightCoord{"height coordinate", this};

  /// Ensure calculated DateTimes lie within the observation window.
  oops::Parameter<bool> keep_in_window{"keep in window", false, this};

  /// Require pressure to be sorted in descending order.
  /// Only set this to `false` if you are happy for the data to potentially
  /// be unsorted, which could lead to unexpected behaviour.
  oops::Parameter<bool> RequireDescendingPressureSort
    {"require descending pressure sort", true, this};
};

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
  typedef Cal_ProfileHorizontalDriftParameters Parameters_;

  Cal_ProfileHorizontalDrift(const Parameters_ &options,
                             const ObsFilterData &data,
                             const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                             const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  /// Height coordinate name.
  std::string heightCoord_;

  /// Keep observations in window:
  bool keep_in_window_;

  /// Require pressure to be sorted in descending order.
  bool requireDescendingPressureSort_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_PROFILEHORIZONTALDRIFT_H_
