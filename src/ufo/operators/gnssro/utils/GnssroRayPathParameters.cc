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
#include <cmath>
#include <memory>
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"

namespace ufo {

static constexpr double KM_TO_METERS = 1000.0;

// -----------------------------------------------------------------------------
GnssroRayPathParameters::GnssroRayPathParameters()
  : rayPathGenType_(DEFAULT_RAY_PATH_GEN_TYPE)
  , approxRayLength_(DEFAULT_APPROX_RAY_LENGTH_KM * KM_TO_METERS)
  , expectedRayLength_(0.0)
  , top2D_(DEFAULT_TOP2D_KM * KM_TO_METERS)
  , horizontalRes_(DEFAULT_HORIZONTAL_RES_KM * KM_TO_METERS)
  , nHoriz_(DEFAULT_NHORIZ)
{
  computeDerivedParameters();
}
// -----------------------------------------------------------------------------

GnssroRayPathParameters::GnssroRayPathParameters(
    const std::string & rayPathGenType,
    double approxRayLengthKm,
    double top2DKm,
    double resKm,
    int nHoriz
) : rayPathGenType_(rayPathGenType)
  , approxRayLength_(approxRayLengthKm * KM_TO_METERS)
  , expectedRayLength_(0.0)
  , top2D_(top2DKm * KM_TO_METERS)
  , horizontalRes_(resKm * KM_TO_METERS)
  , nHoriz_(nHoriz)
{
  computeDerivedParameters();
}
// -----------------------------------------------------------------------------

void GnssroRayPathParameters::computeDerivedParameters()
{
  if (approxRayLength_ > 0.0) {
    // Compute nHoriz_ if max_ray_length was specified.
    int originalNHoriz = nHoriz_;

    // Ensure nHoriz_ is an odd number and within reasonable bounds.
    int nHorizFloor = static_cast<int>(std::floor(approxRayLength_ / horizontalRes_));
    int nHorizCeil  = static_cast<int>(std::ceil(approxRayLength_ / horizontalRes_));
    if (nHorizFloor % 2 != 0) {        // Floor is odd
      nHoriz_ = nHorizFloor;
    } else if (nHorizCeil % 2 != 0) {  // Ceiling is odd
      nHoriz_ = nHorizCeil;
    } else {  // Both ceiling and floor are even (original number was even int)
      nHoriz_ = nHorizFloor + 1;
    }

    // Emit warning if we changed the configured n_horiz.
    if (originalNHoriz != DEFAULT_NHORIZ && nHoriz_ != DEFAULT_NHORIZ) {
      oops::Log::warning() << "GnssroRayPathParameters: overriding configured nHoriz="
              << originalNHoriz << " with computed value=" << nHoriz_
              << " from approxRayLength=" << approxRayLength_
              << " and horizontalRes=" << horizontalRes_ << std::endl;
    }

    // Ensure nHoriz_ is an odd number and within reasonable bounds.
    if (nHoriz_ <= 0) {
      oops::Log::error() << "GnssroRayPathParameters: adjusting nHoriz="
              << nHoriz_ << " up to 1; calculated from approxRayLength="
              << approxRayLength_ << " and horizontalRes=" << horizontalRes_
              << std::endl;
      nHoriz_ = 1;
    }
    if (nHoriz_ > MAX_NHORIZ) {
      oops::Log::error() << "GnssroRayPathParameters: adjusting nHoriz="
              << nHoriz_ << " down to " << MAX_NHORIZ
              << "; calculated from approxRayLength="
              << approxRayLength_ << " and horizontalRes=" << horizontalRes_
              << std::endl;
      nHoriz_ = MAX_NHORIZ;
    }

    oops::Log::debug() << "GnssroRayPathParameters: nHoriz=" << nHoriz_
            << " from approxRayLength=" << approxRayLength_ << " and horizontalRes="
            << horizontalRes_ << std::endl;
  } else {
    // ray length not configured; use nHoriz
    // Ensure nHoriz_ is an odd number and within reasonable bounds.
    if (nHoriz_ <= 0 || nHoriz_ > MAX_NHORIZ || nHoriz_ % 2 == 0) {
      throw eckit::BadValue(
          "GnssroRayPathParameters: Bad value for nhoriz: "
          + std::to_string(nHoriz_)
          + "; must be an odd number between 1 and "
          + std::to_string(MAX_NHORIZ) + " (inclusive)", Here());
    }
  }

  // Compute expected ray length
  expectedRayLength_ = nHoriz_ * horizontalRes_;
  return;
}
// -----------------------------------------------------------------------------

}  // namespace ufo

