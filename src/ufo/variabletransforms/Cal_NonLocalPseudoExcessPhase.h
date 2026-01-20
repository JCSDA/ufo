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
#ifndef UFO_VARIABLETRANSFORMS_CAL_NONLOCALPSEUDOEXCESSPHASE_H_
#define UFO_VARIABLETRANSFORMS_CAL_NONLOCALPSEUDOEXCESSPHASE_H_

#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ufo/operators/gnssro/utils/GnssroProfileRayPaths.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathOrchestrator.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the wind components transform.
class Cal_NonLocalPseudoExcessPhaseParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_NonLocalPseudoExcessPhaseParameters,
                           VariableTransformParametersBase);

 public:
  /// Observation group name. Default is ObsValue.
  oops::Parameter<std::string> group{"group", "ObsValue", this};
  oops::Parameter<std::string> refrVariable{"refractivity variable",
                                            "atmosphericRefractivity", this};
  oops::Parameter<std::string> nlpepVariable{"nlpep variable",
                                             "nonLocalPseudoExcessPhase", this};
  oops::Parameter<std::string> rayPathGenType{"ray_path_gen_type",
          GnssroRayPathParameters::DEFAULT_RAY_PATH_GEN_TYPE, this};
  oops::Parameter<double> approxRayLength{"ray_length",
          GnssroRayPathParameters::DEFAULT_APPROX_RAY_LENGTH_KM, this};
  oops::Parameter<double> res{"res",
          GnssroRayPathParameters::DEFAULT_HORIZONTAL_RES_KM, this};
  oops::Parameter<double> top2D{"top_2d",
          GnssroRayPathParameters::DEFAULT_TOP2D_KM, this};
  oops::Parameter<int> nHoriz{"n_horiz",
          GnssroRayPathParameters::DEFAULT_NHORIZ, this};
};

/*!
* \brief nonLocalPseudoExcessPhase filter
*
* \details  Performs a variable conversion from atmosphericRefractivity to nonLocalPseudoExcessPhase.
*  nonLocalPseudoExcessPhase is a path-integrated form of refractivity, dependent on
*  a pre-defined ray path. This ray path must match the ray path used by an associated
*  GNSSRO forward operator, e.g. RefNLPEP2D.
*  The newly calculated variables are included in the same obs space.
///
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_NonLocalPseudoExcessPhase: public TransformBase {
 public:
  typedef Cal_NonLocalPseudoExcessPhaseParameters Parameters_;

  Cal_NonLocalPseudoExcessPhase(const Parameters_ &options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  // Member data
  std::string refrvariable_;
  std::string nlpepvariable_;
  std::string group_;
  GnssroRayPathParameters rpParams_;
  GnssroRayPathOrchestrator orchestrator_;
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_NONLOCALPSEUDOEXCESSPHASE_H_
