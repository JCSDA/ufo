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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATH_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATH_H_

#include <cassert>
#include <cstddef>  // For std::size_t
#include <cstdint>  // For uint32_t
#include <vector>
#include "ufo/operators/gnssro/utils/GnssroRayTrajectory.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Class to represent represent the path of a single GNSSRO ray.
/// Typically, each sample in a set of RO observations becomes a ray.
class GnssroRayPath
{
 public:
  typedef std::vector<double> DoubleArray;
  typedef std::vector<float> RealArray;
  typedef std::vector<std::uint32_t> UIntArray;

  GnssroRayPath(uint32_t num_nodes, float tpLat, float tpLon, float azimuth,
      float tpAlt, bool assumeSphericalSymmetry = false);
  GnssroRayPath() = default;
  GnssroRayPath(const GnssroRayPath & other) = delete;
  GnssroRayPath(GnssroRayPath && other) noexcept = delete;
  GnssroRayPath & operator=(const GnssroRayPath & other) = delete;
  GnssroRayPath & operator=(GnssroRayPath && other) noexcept = delete;
  ~GnssroRayPath();

  // Accessors.
  std::uint32_t centerNodeIdx() const {return segmentLengths_.size() / 2;}
  std::uint32_t numEdges() const {return geomHeights_.size();}
  std::uint32_t numNodes() const {return segmentLengths_.size();}
  std::uint32_t numLocations() const {
    return (assumeSphericalSymmetry_ ? 1 : numNodes());
  }
  std::uint32_t edgesInHalfRay() const {return numEdges() / 2;}
  bool assumeSphericalSymmetry() const {return assumeSphericalSymmetry_;}
  const GnssroRayTrajectory& traj() const {return traj_;}
  GnssroRayTrajectory& traj() {return traj_;}

  // Read-only accessors to array data.
  float geomHgt(uint32_t edgeIdx) const {return geomHeights_[edgeIdx];}
  float lat(uint32_t nodeIdx) const {
    return (assumeSphericalSymmetry_ ? tpLat_ : latitudes_[nodeIdx]);
  }
  float lon(uint32_t nodeIdx) const {
    return (assumeSphericalSymmetry_ ? tpLon_ : longitudes_[nodeIdx]);
  }
  double segLen(uint32_t nodeIdx) const {return segmentLengths_[nodeIdx];}

  // Read only accessors for the vectors
  const RealArray& lats() const {return latitudes_;}
  const RealArray& lons() const {return longitudes_;}

  // Read-write accessors to array data.
  float & geomHgt(uint32_t edgeIdx) {
    assert(edgeIdx < geomHeights_.size());
    return geomHeights_[edgeIdx];
  }
  float & lat(uint32_t nodeIdx) {
    if (assumeSphericalSymmetry_) {
      return tpLat_;
    } else {
      assert(nodeIdx < latitudes_.size());
      return latitudes_[nodeIdx];
    }
  }
  float & lon(uint32_t nodeIdx) {
    if (assumeSphericalSymmetry_) {
      return tpLon_;
    } else {
      assert(nodeIdx < longitudes_.size());
      return longitudes_[nodeIdx];
    }
  }
  double & segLen(uint32_t nodeIdx) {
    assert(nodeIdx < segmentLengths_.size());
    return segmentLengths_[nodeIdx];
  }

  // Remove elements from front and back of vector so there are edgeOffset points
  // around the center.
  void removeUnusedNodes(std::uint32_t edgeOffset);

  // Resize to specified number of nodes by removing elements from the end,
  // if necessary.
  void shrinkNodesTo(std::uint32_t newNumNodes);

  // Method to compute the total length from the segment lengths.
  void updateTotalLength();

  float tpLat() const {return tpLat_;}
  float tpLon() const {return tpLon_;}
  float azimuth() const {return azimuth_;}
  float tpAlt() const {return tpAlt_;}
  float totalLength() const {return totalLength_;}
  std::uint32_t tpNodeIdx() const {return tpNodeIdx_;}
  void setTpNodeIdx(std::uint32_t tpNodeIdx) {tpNodeIdx_ = tpNodeIdx;}
  std::uint32_t beginLocIdx() const {return beginLocIdx_;}
  std::uint32_t endLocIdx() const {return endLocIdx_;}
  void setLocIdxRange(std::uint32_t beginLocIdx, std::uint32_t endLocIdx);

 private:
  //  Member data.
  float tpLat_;                   // Latitude of tangent point defining this ray.
  float tpLon_;                   // Longitude of tangent point defining this ray.
  float azimuth_;                 // Azimuth angle defining this ray.
  float tpAlt_;                   // Altitude of tangent point (Geometric height above MSL in m)
  double totalLength_;            // sum of all segment lengths.
  uint32_t tpNodeIdx_;            // Node index of the tangent point.
  bool assumeSphericalSymmetry_;  // Each node with use the tangent point lat, lon
  RealArray geomHeights_;         // For bottom edges of nodes (m). Last edge is top of highest node
  RealArray latitudes_;           // For edges (degrees)
  RealArray longitudes_;          // For edges (degrees)
  DoubleArray segmentLengths_;    // Length of segments centered on node (m)
  std::uint32_t beginLocIdx_;     // begin index for pathsGroupedByLocation
  std::uint32_t endLocIdx_;       // end index for pathsGroupedByLocation
  GnssroRayTrajectory traj_;      // 2D Trajectory used in TLAD
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATH_H_
