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
#include <algorithm>  // for std::copy, std::nth_element
#include <cassert>
#include <limits>     // for NaN
#include <utility>
#include "oops/util/Logger.h"
#include "ufo/operators/gnssro/utils/GnssroProfileRayPaths.h"
#include "ufo/operators/gnssro/utils/GnssroRayPath.h"

namespace ufo {

// Helper function to approximate a median calculation.
// This sorts the input vector just enough to determine the central value.
// However, for vectors with an even number of elements, it does
// not take the mean of the two values closest to the center.
static float approx_median(
    const std::vector<float> & v,
    std::size_t start_offset,
    std::size_t end_offset,
    std::vector<float> & tmp_v)
{
  // Copy the subset of the source vector to the temporary vector.
  std::copy(v.begin() + start_offset, v.begin() + end_offset, tmp_v.begin());

  // Determine approximate location of the median in sorted order.
  // If the size is even, the median should be the average of the
  // two center values. However, we don't need to be that precise.
  std::size_t median_loc = tmp_v.size() / 2;

  // Use nth_element to sort the data just enough to place
  // the value for the median location in the right place.
  // This is more efficient than sorting the entire vector.
  std::nth_element(tmp_v.begin(), tmp_v.begin() + median_loc, tmp_v.end());
  return tmp_v[median_loc];
}

// -----------------------------------------------------------------------------
/// Class to represent represent the path of a single GNSSRO ray.
/// Typically, each sample in a set of RO observations becomes a ray.
GnssroProfileRayPaths::GnssroProfileRayPaths()
  : seqNum_(-1)
  , profileLat_(std::numeric_limits<float>::quiet_NaN())
  , profileLon_(std::numeric_limits<float>::quiet_NaN())
  , profileAzimuth_(std::numeric_limits<float>::quiet_NaN())
  , rays_()
{
  oops::Log::trace() << "GnssroProfileRayPaths default c'tor" << std::endl;
}
// -----------------------------------------------------------------------------

GnssroProfileRayPaths::GnssroProfileRayPaths(
    const GnssroProfileSlice & slice,
    const std::vector<float> & lats,
    const std::vector<float> & lons,
    const std::vector<float> & azimuths)
  : seqNum_(slice.seqNum())
  , profileLat_(std::numeric_limits<float>::quiet_NaN())
  , profileLon_(std::numeric_limits<float>::quiet_NaN())
  , profileAzimuth_(std::numeric_limits<float>::quiet_NaN())
  , rays_()
{
  oops::Log::trace() << "GnssroProfileRayPaths parameterized c'tor" << std::endl;
  setProfileValues(slice, lats, lons, azimuths);

  // The size of the profile is the maximum number of rays we will have,
  // assuming no observations are filtered out.
  rays_.reserve(slice.size());
}
// -----------------------------------------------------------------------------

GnssroProfileRayPaths::~GnssroProfileRayPaths()
{
}
// -----------------------------------------------------------------------------

void GnssroProfileRayPaths::setProfileValues(
    const GnssroProfileSlice & slice,
    const std::vector<float> & lats,
    const std::vector<float> & lons,
    const std::vector<float> & azimuths)
{
  oops::Log::trace() << "GnssroProfileRayPaths::setProfileValues - enter" << std::endl;
  if (lats.size() == 0) {
      return;
  }

  // Define a scratch vector, then compute an approximate median for each array.
  std::vector<float> tmp_v(slice.size());
  profileLat_ = approx_median(lats, slice.start(), slice.end(), tmp_v);
  profileLon_ = approx_median(lons, slice.start(), slice.end(), tmp_v);
  profileAzimuth_ = approx_median(azimuths, slice.start(), slice.end(), tmp_v);

  // Lats, lons, and azimuths should all have the same size (number of samples
  // in all the obs passed to this rank).  The slice size is the number of samples
  // in a single profile, which must be <= this size.
  assert(lons.size() == lats.size());
  assert(azimuths.size() == lats.size());
  assert(slice.size() <= lats.size());
  oops::Log::debug() << "GnssroProfileRayPaths seqNum=" << slice.seqNum() << ": lat "
      << *(std::min_element(lats.cbegin() + slice.start(), lats.cbegin() + slice.end())) << " to "
      << *(std::max_element(lats.cbegin() + slice.start(), lats.cbegin() + slice.end())) << " ("
      << profileLat_ << "), lon "
      << *(std::min_element(lons.cbegin() + slice.start(), lons.cbegin() + slice.end())) << " to "
      << *(std::max_element(lons.cbegin() + slice.start(), lons.cbegin() + slice.end())) << " ("
      << profileLon_ << "), azimuth "
      << *(std::min_element(azimuths.cbegin() + slice.start(), azimuths.cbegin() + slice.end()))
      << " to "
      << *(std::max_element(azimuths.cbegin() + slice.start(), azimuths.cbegin() + slice.end()))
      << " (" << profileAzimuth_ << ")" << std::endl;
  oops::Log::trace() << "GnssroProfileRayPaths::setProfileValues - exit" << std::endl;
}
// -----------------------------------------------------------------------------

void GnssroProfileRayPaths::addRay(std::unique_ptr<GnssroRayPath> && ray) noexcept
{
  rays_.push_back(std::move(ray));
}
// -----------------------------------------------------------------------------

// Count all nodes in all rays in the profile.
std::size_t GnssroProfileRayPaths::totalNodes() const
{
  std::size_t nodeCnt = 0;
  for (const_iterator rayItr = cbegin(); rayItr != cend(); ++rayItr)
  {
    const GnssroProfileRayPaths::Ray_& ray = *rayItr;
    nodeCnt += ray->numNodes();
  }
  return nodeCnt;
}
// -----------------------------------------------------------------------------

// Count all sampled locations in all rays in the profile.
std::size_t GnssroProfileRayPaths::totalSampledLocations() const
{
  std::size_t sampLocs = 0;
  for (const_iterator rayItr = cbegin(); rayItr != cend(); ++rayItr)
  {
    const GnssroProfileRayPaths::Ray_& ray = *rayItr;
    sampLocs += ray->numLocations();
  }
  return sampLocs;
}
// -----------------------------------------------------------------------------


}  // namespace ufo

