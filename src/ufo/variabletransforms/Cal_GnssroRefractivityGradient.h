/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#ifndef UFO_VARIABLETRANSFORMS_CAL_GNSSROREFRACTIVITYGRADIENT_H_
#define UFO_VARIABLETRANSFORMS_CAL_GNSSROREFRACTIVITYGRADIENT_H_

#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ufo/operators/gnssro/utils/GnssroProfileExtractor.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the wind components transform.
class Cal_GnssroRefractivityGradientParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_GnssroRefractivityGradientParameters,
                           VariableTransformParametersBase);

 public:
  /// Observation group name. Default is ObsValue.
  oops::Parameter<std::string> group{"group", "ObsValue", this};
  oops::Parameter<std::string> refrVariable{"refractivity variable",
                                            "atmosphericRefractivity", this};
  oops::Parameter<std::string> heightVariable{"height variable",
                                            "height", this};
  oops::Parameter<std::string> gradVariable{"gradient variable",
                                            "atmosphericRefractivityGradient", this};
  oops::Parameter<float> minSuperRefraction{"min super refraction", 0.5, this};
  oops::Parameter<bool> calculateDuctingFlag{"calculate ducting flag", false, this};
  oops::Parameter<bool> calculateProfileDuctingFlag{"calculate profile ducting flag",
                                                    false, this};
};

/*!
* \brief GNSSRO refractivity gradient filter
*
* \details  Performs the vertical derivative of refractivity wrt geometric height dN/dz.
*  The newly calculated variables are included in the same obs space in a group with
*  "Derived" as the prefix; default is DerivedObsValue.
///
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_GnssroRefractivityGradient: public TransformBase {
 public:
  static constexpr int NONDUCTING = 0;
  static constexpr int MIN_SUPER_REFRACTION = 1;
  static constexpr int NEAR_DUCTING = 2;
  static constexpr int QUITE_NEAR_DUCTING = 3;
  static constexpr int VERY_NEAR_DUCTING = 4;
  static constexpr int EXTREMELY_NEAR_DUCTING = 5;
  static constexpr int DUCTING = 6;

  typedef Cal_GnssroRefractivityGradientParameters Parameters_;

  Cal_GnssroRefractivityGradient(const Parameters_ &options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  // Member data
  std::string refrvariable_;
  std::string heightvariable_;
  std::string gradvariable_;
  std::string group_;
  float minSuperRefraction_;    // Smallest fraction of critical value tracked.
  float minSRThreshold_;        // units: N-units / meter
  float nearThreshold_;
  float quiteNearThreshold_;
  float veryNearThreshold_;
  float extremelyNearThreshold_;
  bool calcDuctingFlag_;
  bool calcProfileDuctingFlag_;
  GnssroProfileExtractor extractor_;
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_GNSSROREFRACTIVITYGRADIENT_H_
