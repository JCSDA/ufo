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
#ifndef UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DPARAMETERS_H_
#define UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DPARAMETERS_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the GnssroRefNLPEP2D operator.
class ObsGnssroRefNLPEP2DOptions : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGnssroRefNLPEP2DOptions, Parameters)

 public:
  /// \brief GnssroRayPathGenerator to use to create ray path.
  oops::Parameter<std::string> rayPathGenType{"ray_path_gen_type",
          GnssroRayPathParameters::DEFAULT_RAY_PATH_GEN_TYPE, this};

  /// \brief Approximate ray length in km.
  /// Ray length and res will be used to override the n_horiz value if ray_length
  /// is set to a value > 0. The effective ray length may be adjusted slightly lower
  /// or higher to meet the requirement of some ray path generators, such as having
  /// an odd number of nodes.
  oops::Parameter<double> approxRayLength{"ray_length",
          GnssroRayPathParameters::DEFAULT_APPROX_RAY_LENGTH_KM, this};

  /// \brief horizontal resolution of nodes in ray path in km.
  /// This should be close to the native resolution of the NWP model.
  oops::Parameter<double> res{"res",
          GnssroRayPathParameters::DEFAULT_HORIZONTAL_RES_KM, this};

  /// \brief The highest geometric height in km that will use 2D rays.
  /// Levels above this height will be treated as 1D rays (one node at tangent point,
  /// assuming the atmosphere is spherically symmetric around tangent point).
  oops::Parameter<double> top2D{"top_2d",
          GnssroRayPathParameters::DEFAULT_TOP2D_KM, this};

  /// \brief number of nodes in each 2D ray. Must be odd.
  oops::Parameter<int> nHoriz{"n_horiz",
          GnssroRayPathParameters::DEFAULT_NHORIZ, this};

  /// \brief The algorithm to use for computing refractivity and its derivatives.
  oops::Parameter<std::string> refrAlgorithm{"refr_algo", "RuegerBevis", this};

  /// \brief true if the refractivity algorithm should assume a compressible atmosphere.
  /// Set to false for an uncompressible atmosphere.
  oops::Parameter<bool> useCompress{"use_compress", true, this};
};

class ObsGnssroRefNLPEP2DParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsGnssroRefNLPEP2DParameters, ObsOperatorParametersBase)
 public:
  oops::Parameter<ObsGnssroRefNLPEP2DOptions> options{"obs options", {}, this};
};


}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DPARAMETERS_H_
