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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILERAYPATHS_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILERAYPATHS_H_

#include <cstddef>  // For std::size_t
#include <memory>
#include <vector>
#include "ufo/operators/gnssro/utils/GnssroProfileSlice.h"
#include "ufo/operators/gnssro/utils/GnssroRayPath.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Class to represent represent the ray paths associated with a single GNSSRO
/// profile.
class GnssroProfileRayPaths
{
 public:
  typedef std::unique_ptr<GnssroRayPath> Ray_;
  typedef std::vector<Ray_> Rays_;
  typedef std::vector<Ray_>::const_iterator const_iterator;
  typedef std::vector<Ray_>::const_reverse_iterator const_reverse_iterator;

  GnssroProfileRayPaths();
  GnssroProfileRayPaths(const GnssroProfileSlice & slice,
                        const std::vector<float> & lats,
                        const std::vector<float> & lons,
                        const std::vector<float> & azimuths);
  ~GnssroProfileRayPaths();
  GnssroProfileRayPaths(const GnssroProfileRayPaths & other) = delete;
  GnssroProfileRayPaths(GnssroProfileRayPaths && other) noexcept = delete;
  GnssroProfileRayPaths & operator=(const GnssroProfileRayPaths & other) = delete;
  GnssroProfileRayPaths & operator=(GnssroProfileRayPaths && other) noexcept = delete;

  void setProfileValues(const GnssroProfileSlice & slice,
                        const std::vector<float> & lats,
                        const std::vector<float> & lons,
                        const std::vector<float> & azimuths);

  //  Ray accessors.
  void addRay(std::unique_ptr<GnssroRayPath> && ray) noexcept;
  std::size_t numRays() const {return rays_.size();}
  const Ray_ & ray(std::size_t rayIdx) const {return rays_[rayIdx];}
  Ray_ & ray(std::size_t rayIdx) {return rays_[rayIdx];}

  // Ray iterators
  const Rays_& rays() const {return rays_;}
  const_iterator cbegin() const {return rays_.cbegin();}
  const_iterator cend() const {return rays_.cend();}
  const_reverse_iterator crbegin() const {return rays_.crbegin();}
  const_reverse_iterator crend() const {return rays_.crend();}

  //  Profile metadata accessors.
  int   seqNum() const {return seqNum_;}
  float profileLat() const {return profileLat_;}
  float profileLon() const {return profileLon_;}
  float profileAzimuth() const {return profileAzimuth_;}

  //  Aggregate quantities over all rays in RO profile
  std::size_t totalNodes() const;
  std::size_t totalSampledLocations() const;

 private:
  // These are representative values for the entire profile (all in degrees).
  int   seqNum_;          // Profile identifier: unique within an input obs file.
  float profileLat_;
  float profileLon_;
  float profileAzimuth_;
  Rays_ rays_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILERAYPATHS_H_
