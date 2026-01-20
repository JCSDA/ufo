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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHPARAMETERS_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHPARAMETERS_H_

#include <string>

namespace ufo {

// -----------------------------------------------------------------------------
class GnssroRayPathParameters {
 public:
  //  Class to hold configurable parameters for GnssroRayPathGenerator and
  //  GnssroRayPathOrchestrator classes.
  //  This should hold the super set of parameters needed by all subclasses.
  static constexpr const char DEFAULT_RAY_PATH_GEN_TYPE[] = "StraightLine";
  static constexpr double DEFAULT_APPROX_RAY_LENGTH_KM = -1.0;
  static constexpr double DEFAULT_TOP2D_KM = 90.0;
  static constexpr double DEFAULT_HORIZONTAL_RES_KM = 11.0;
  static constexpr int DEFAULT_NHORIZ = 11;
  static constexpr int MAX_NHORIZ = 67;

  GnssroRayPathParameters();
  GnssroRayPathParameters(
      const std::string & rayPathGenType,
      double approxRayLengthKm,
      double top2DKm,
      double resKm,
      int nHoriz);
  ~GnssroRayPathParameters() {}
  GnssroRayPathParameters(const GnssroRayPathParameters & other) = default;
  GnssroRayPathParameters & operator=(const GnssroRayPathParameters & other) = default;

  void computeDerivedParameters();

  const std::string & rayPathGenType() const {return rayPathGenType_;}
  double approxRayLength() const {return approxRayLength_;}
  double expectedRayLength() const {return expectedRayLength_;}
  double top2D() const {return top2D_;}
  double horizontalRes() const {return horizontalRes_;}
  int nHoriz() const {return nHoriz_;}

 private:
  //  Member data
  std::string rayPathGenType_;  // Name for GnssroRayPathGenerator subclass
  double approxRayLength_;      // Configured ray length (m)
  double expectedRayLength_;    // ray length (m) resolved from n_horiz and res
  double top2D_;                // Highest geom height (m) for 2D rays
  double horizontalRes_;        // Horizontal resolution (m)
  int nHoriz_;                  // Number of nodes in each 2D ray
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHPARAMETERS_H_
