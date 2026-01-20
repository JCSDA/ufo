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
// Define compilation symbol DEBUG_GNSSRO_RAY_PATH (e.g. with -D compiler option)
// to enable debug messages in this file.
#include <algorithm>
#include <limits>
#include <numeric>
#include "oops/util/Logger.h"
#include "ufo/operators/gnssro/utils/GnssroRayPath.h"

namespace ufo {

template <typename T>
static void eraseVectorFromBothEnds(std::vector<T> & v, std::uint32_t numNodesToRemove)
{
  // Erase from end first, since this does not require moving values.
  v.erase(v.end() - numNodesToRemove, v.end());
  // Erase from beginning second.
  v.erase(v.begin(), v.begin() + numNodesToRemove);
}

// -----------------------------------------------------------------------------
/// Class to represent represent the path of a single GNSSRO ray.
/// Typically, each sample in a set of RO observations becomes a ray.
GnssroRayPath::GnssroRayPath(
    std::uint32_t numNodes,
    float tpLat,
    float tpLon,
    float azimuth,
    float tpAlt,
    bool assumeSphericalSymmetry)
  : tpLat_(tpLat)
  , tpLon_(tpLon)
  , azimuth_(azimuth)
  , tpAlt_(tpAlt)
  , totalLength_(0.0)
  , tpNodeIdx_(std::numeric_limits<std::uint32_t>::max())
  , assumeSphericalSymmetry_(assumeSphericalSymmetry)
  , geomHeights_(numNodes + 1)  // Dimensioned by edges
  , latitudes_()                // Dimensioned by nodes (not used w/ SphericalSym)
  , longitudes_()               // Dimensioned by nodes (not used w/ SphericalSym)
  , segmentLengths_(numNodes)   // Dimensioned by nodes (center points)
  , beginLocIdx_(0)
  , endLocIdx_(0)
  , traj_()
{
  // Only allocate lat-lon vectors if we are not assuming spherical symmetry.
  // Otherwise, we assume all nodes are located at the tangent point.
  if (!assumeSphericalSymmetry_) {
    latitudes_.assign(numNodes, 0.0);
    longitudes_.assign(numNodes, 0.0);
  }
}

// -----------------------------------------------------------------------------

GnssroRayPath::~GnssroRayPath()
{
}
// -----------------------------------------------------------------------------

void GnssroRayPath::removeUnusedNodes(std::uint32_t edgeOffset)
{
  // Compute nodes to remove from the front and the back of each vector.
  std::uint32_t tpIdx = centerNodeIdx();
  std::uint32_t numNodesToRemove = tpIdx - edgeOffset;

  if (numNodesToRemove == 0)
    return;

  eraseVectorFromBothEnds(geomHeights_, numNodesToRemove);
  eraseVectorFromBothEnds(segmentLengths_, numNodesToRemove);
  if (!assumeSphericalSymmetry_) {
    eraseVectorFromBothEnds(latitudes_, numNodesToRemove);
    eraseVectorFromBothEnds(longitudes_, numNodesToRemove);

    // Check for off-by-one error loading these arrays.
    if (geomHeights_.size() > 0 && latitudes_[0] == 0.0 && longitudes_[0] == 0.0
        && geomHeights_[0] == 0.0)
    {
      oops::Log::error() << "GnssroRayPath: first lat/lon/hgt == 0.0 after erasing numNodes="
          << numNodesToRemove << " from each side of tp to shrink ray to " << latitudes_.size()
          << " nodes" << std::endl;
    }
  }

  // Reset the tangent point node index. Since we removed nodes symmetrically
  // from both ends, the tangent point will still be the center node.
  setTpNodeIdx(centerNodeIdx());
}
// -----------------------------------------------------------------------------

void GnssroRayPath::shrinkNodesTo(std::uint32_t newNumNodes)
{
  if (newNumNodes < numNodes()) {
    uint32_t numToRemove = numNodes() - newNumNodes;
    geomHeights_.erase(geomHeights_.end() - numToRemove, geomHeights_.end());
    segmentLengths_.erase(segmentLengths_.end() - numToRemove, segmentLengths_.end());
    if (!assumeSphericalSymmetry_) {
      latitudes_.erase(latitudes_.end() - numToRemove, latitudes_.end());
      longitudes_.erase(longitudes_.end() - numToRemove, longitudes_.end());
    }
  } else if (newNumNodes > numNodes()) {
    oops::Log::error() << "GnssroRayPath::shrinkNodesTo: newNumNodes=" << newNumNodes
            << " > oldNumNodes=" << numNodes() << "; not changing ray size" << std::endl;
  }
  // Compute the total length of the ray. We only removed nodes after the last used node,
  // so the tangent point node will not have changed.
  updateTotalLength();
  return;
}
// -----------------------------------------------------------------------------

void GnssroRayPath::updateTotalLength()
{
  double newLength = std::accumulate(segmentLengths_.begin(), segmentLengths_.end(), 0.0);
  totalLength_ = newLength;
}
// -----------------------------------------------------------------------------

void GnssroRayPath::setLocIdxRange(std::uint32_t beginLocIdx, std::uint32_t endLocIdx)
{
  beginLocIdx_ = beginLocIdx;
  endLocIdx_ = endLocIdx;
}
// -----------------------------------------------------------------------------

}  // namespace ufo

