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
#include <cassert>
#include <cmath>
#include <cstddef>  // for std::size_t
#include <cstdlib>
#include <memory>
#include <string>
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/operators/gnssro/utils/GnssroRayPath.h"
#include "ufo/operators/gnssro/utils/StraightLineRayPathGenerator.h"

namespace ufo {

static const double RAD2DEG = 180.0 / M_PI;
static const double DEG2RAD = M_PI / 180.0;

template <typename T>
static T normalizeLongitude_180(T lon)
{
  // Normalize longitude in degrees to -180 to 180 range.
  T angle = std::fmod(lon + static_cast<T>(180), static_cast<T>(360));
  if (angle < static_cast<T>(0))
    angle += static_cast<T>(360);
  angle -= static_cast<T>(180);
  return angle;
}

template <typename T>
static T normalizeLongitude_360(T lon)
{
  static const T three_sixty = static_cast<T>(360);
  return std::fmod(std::fmod(lon, three_sixty) + three_sixty, three_sixty);
}

/*
 * Geometry for straightline calculations.
 *                  d = rayDist
 *                -----
 *                \   |
 *  edgeRadius = R \  | R0 = radius0
 *                  \ |
 *                   \|
 *  Angle between R and R0 is baseTheta.
 */
StraightLineRayPathGenerator::NodeCalculator::NodeCalculator(
        double tpGeomHeight,
        double localEarthRadius,
        float azimuth,
        float tpLat,
        float tpLon,
        const GnssroRayPathParameters & params)
  : tpGeomHeight_(tpGeomHeight)
  , localEarthRadius_(localEarthRadius)
  , radius0_(tpGeomHeight_ + localEarthRadius_)
  , radius0_sq_(radius0_ * radius0_)
  , azimuth_((azimuth > 180.0) ? azimuth - 360.0 : azimuth)
  , tpLat_(tpLat)
  , tpLon_(tpLon)
  , params_(params)
  , azimuthRad_(azimuth_ * DEG2RAD)
  , tpLatRad_(tpLat_ * DEG2RAD)
  , tpLonRad_(tpLon_ * DEG2RAD)
  , cosAzimuth_(std::cos(azimuthRad_))
  , signAzimuth_((azimuth_ >= 0.0) ? 1 : -1)
  , cosTpLat_(std::cos(tpLatRad_))
  , cosTpLon_(std::cos(tpLonRad_))
  , sinTpLat_(std::sin(tpLatRad_))
  , sinTpLon_(std::sin(tpLonRad_))
  , edgeIncrement_(0.5 * params_.horizontalRes())  // increment for rayDist
  , edgeDist_(0.0)                                // Length of half ray (from TP)
  , edgeRadius_(radius0_)                         // Radius to current edge
  , baseTheta_(0.0)                               // Angle to current node (not edge)
{
}

StraightLineRayPathGenerator::NodeCalculator::~NodeCalculator()
{
}

// Compute increment ray distance and recompute quantities that depend on it.
void StraightLineRayPathGenerator::NodeCalculator::extendRay()
{
  // Compute the distance to the next edge between nodes.
  edgeDist_ += edgeIncrement_;

  // After the first ray distance computation, we increment by full resolution.
  edgeIncrement_ = params_.horizontalRes();

  // Compute the edge radius, used to determine geometric height of an edge
  edgeRadius_ = std::sqrt(edgeDist_ * edgeDist_ + radius0_sq_);

  // Compute the angle to the current node, which is half of a segment length
  // closer to the tangent point than the current edge.
  // This angle is 0 at the tangent point.
  double nodeDist = edgeDist_ - 0.5 * edgeIncrement_;
  baseTheta_ = std::atan2(nodeDist, radius0_);
  return;
}

std::pair<float, float> StraightLineRayPathGenerator::NodeCalculator::lonLatAlongRay(double theta)
{
  // This algorithm is a copy of the Fortran logic from ROPP, except we use a
  // variable theta, rather than a constant value based on distances at the equator.
  double sinTheta = std::sin(theta);
  double cosTheta = std::cos(theta);

  double sinLat = sinTheta * cosAzimuth_ * cosTpLat_ + sinTpLat_ * cosTheta;
  double latRad = std::asin(sinLat);

  double cosDelta = 1.0;   // For North/South poles
  if (std::abs(latRad) != M_PI / 2.0)  // Not North/South pole
  {
    cosDelta = (cosTheta - (sinTpLat_ * sinLat)) / (cosTpLat_ * std::cos(latRad));
    if (cosDelta > 1.0)
      cosDelta = 1.0;
    if (cosDelta < -1.0)
      cosDelta = -1.0;
  }

  int signTheta = (theta >= 0.0) ? 1 : -1;
  double lonRad;
  if (signTheta == signAzimuth_) {
    lonRad = tpLonRad_ + std::acos(cosDelta);
  } else {
    lonRad = tpLonRad_ - std::acos(cosDelta);
  }

  // Convert lat, lon to degrees; normalize lon to -pi to pi range.
  double latDeg = latRad * RAD2DEG;
  double lonDeg = lonRad * RAD2DEG;
  lonDeg = normalizeLongitude_180(lonDeg);
  return std::pair<double, double>(lonDeg, latDeg);
}

//////////////////////////////////////////////////
// Implementation of StraightLineRayPathGenerator
//////////////////////////////////////////////////
StraightLineRayPathGenerator::StraightLineRayPathGenerator(
    const ioda::ObsSpace & odb,
    const GnssroRayPathParameters & params)
  : GnssroRayPathGenerator(odb, params)
  , obsLat_(odb.nlocs())
  , obsLon_(odb.nlocs())
  , obsAlt_(odb.nlocs())
  , obsRoc_(odb.nlocs())
  , obsUndul_(odb.nlocs())
  , obsAzimuth_(odb.nlocs())
{
  oops::Log::trace() << name() << " entering constructor" << std::endl;

  // Get observation data. This may contain multiple profiles.
  odb.get_db("MetaData", "latitude", obsLat_);
  odb.get_db("MetaData", "longitude", obsLon_);
  odb.get_db("MetaData", "height", obsAlt_);
  odb.get_db("MetaData", "earthRadiusCurvature", obsRoc_);
  odb.get_db("MetaData", "geoidUndulation", obsUndul_);
  odb.get_db("MetaData", "sensorAzimuthAngle", obsAzimuth_);

  oops::Log::trace() << name() << " exiting constructor" << std::endl;
}

StraightLineRayPathGenerator::~StraightLineRayPathGenerator()
{
}

std::unique_ptr<GnssroProfileRayPaths>
StraightLineRayPathGenerator::makeProfileRayPaths(const GnssroProfileSlice & profileSlice)
{
  oops::Log::trace() << name() << "::makeProfileRayPaths - enter" << std::endl;
  auto profileRayPaths = std::make_unique<GnssroProfileRayPaths>(
          profileSlice, obsLat_, obsLon_, obsAzimuth_);

  oops::Log::trace() << name() << "::makeProfileRayPaths - created empty profileRayPaths"
                     << std::endl;
  std::size_t lastObBelowTop = profileSlice.end() - 1;
  std::size_t obIdx = profileSlice.start();
  for (; obIdx <= lastObBelowTop; ++obIdx) {
    // Number of nodes will be constant for all levels.
    // However, lat-lon for all nodes of a ray will be set
    // to the tangent point location for rays above top2D.
    double alt = obsAlt_[obIdx];
    bool assumeSphericalSymmetry = (alt > params_.top2D());
    size_t numNodes = params_.nHoriz();
    oops::Log::debug() << name() << "seqNum=" << profileSlice.seqNum() << ", obIdx=" << obIdx
         << ": checking geometric height=" << alt << " against top2D=" << params_.top2D()
         << ", numNodes=" << numNodes
         << ", assumeSphericalSymmetry=" << assumeSphericalSymmetry
         << std::endl;

    double localEarthRadius = obsRoc_[obIdx] + obsUndul_[obIdx];
    NodeCalculator nodeCalculator(alt,
                                  localEarthRadius,
                                  obsAzimuth_[obIdx],
                                  obsLat_[obIdx],
                                  obsLon_[obIdx],
                                  params_);
    oops::Log::debug() << name() << ": obIdx " << obIdx << ": created NodeCalculator"
         << " with tpLat/Lon " << nodeCalculator.tpLat() << "," << nodeCalculator.tpLon()
         << ", azimuth " << nodeCalculator.azimuth() << std::endl;

    auto rayPath = std::make_unique<GnssroRayPath>(numNodes,
              nodeCalculator.tpLat(),
              nodeCalculator.tpLon(),
              nodeCalculator.azimuth(),
              alt,
              assumeSphericalSymmetry);

    // Set the node index for the tangent point.
    size_t centerNodeIdx = rayPath->centerNodeIdx();
    rayPath->setTpNodeIdx(centerNodeIdx);

    oops::Log::debug() << name() << ": obIdx " << obIdx << ": created GnssroRayPath with "
             << rayPath->numNodes() << " nodes, edgesInHalfRay=" << rayPath->edgesInHalfRay()
             << std::endl;

    for (std::size_t edgeOffset = 0; edgeOffset < rayPath->edgesInHalfRay(); ++edgeOffset)
    {
      nodeCalculator.extendRay();
      double edgeGeomHeight = nodeCalculator.edgeGeomHeight();

      // Set edge-oriented geometric heights, which are symmetric
      // around the tangent point.
      std::size_t rightEdgeIdx = centerNodeIdx + edgeOffset + 1;
      std::size_t leftEdgeIdx = centerNodeIdx - edgeOffset;
      rayPath->geomHgt(rightEdgeIdx) = edgeGeomHeight;
      rayPath->geomHgt(leftEdgeIdx) = edgeGeomHeight;

      //  Set node-oriented values.
      if (edgeOffset == 0) {
        // Tangent point
        assert(rayPath->numNodes() > centerNodeIdx);
        rayPath->segLen(centerNodeIdx) = params_.horizontalRes();
        if (!assumeSphericalSymmetry) {
          rayPath->lat(centerNodeIdx) = nodeCalculator.tpLat();
          rayPath->lon(centerNodeIdx) = nodeCalculator.tpLon();
        }
      } else {
        // Off-tangent points to right and left of tangent point
        std::size_t rightNodeIdx = centerNodeIdx + edgeOffset;
        std::size_t leftNodeIdx = centerNodeIdx - edgeOffset;
        assert(rayPath->numNodes() > rightNodeIdx);
        assert(rayPath->numNodes() > leftNodeIdx);
        rayPath->segLen(rightNodeIdx) = params_.horizontalRes();
        rayPath->segLen(leftNodeIdx) = params_.horizontalRes();

        // Set lat-lons if we are not assuming spherical symmetry.
        if (!assumeSphericalSymmetry) {
          std::pair<float, float> rightLonLat = nodeCalculator.lonLatAlongRay(
                  nodeCalculator.rightTheta());
          std::pair<float, float> leftLonLat = nodeCalculator.lonLatAlongRay(
                  nodeCalculator.leftTheta());
          if (rightLonLat.first == 0.0 && rightLonLat.second == 0.0) {
              oops::Log::warning() << name() << "seqNum=" << profileSlice.seqNum()
                << ", obIdx=" << obIdx << ", edgeIdx=" << rightEdgeIdx << " of "
                << rayPath->numEdges() << ", edgeOffset=" << edgeOffset
                << ": rightLonLat values are both 0.0" << std::endl;
          }
          if (leftLonLat.first == 0.0 && leftLonLat.second == 0.0) {
              oops::Log::warning() << name() << "seqNum=" << profileSlice.seqNum()
                << ", obIdx=" << obIdx << ", edgeIdx=" << leftEdgeIdx << " of "
                << rayPath->numEdges() << ", edgeOffset=" << edgeOffset
                << ": leftLonLat values are both 0.0 " << std::endl;
          }

          rayPath->lon(rightNodeIdx) = rightLonLat.first;
          rayPath->lat(rightNodeIdx) = rightLonLat.second;
          rayPath->lon(leftNodeIdx) = leftLonLat.first;
          rayPath->lat(leftNodeIdx) = leftLonLat.second;
        }  // end of setting lat-lons
      }
    }  // end loop over edges in the ray path.

    rayPath->updateTotalLength();
    oops::Log::debug() << name() << ": obIdx " << obIdx << ": loaded GnssroRayPath with "
            << "numNodes=" << rayPath->numNodes() << ", totalLength=" << rayPath->totalLength()
            << std::endl;
    profileRayPaths->addRay(std::move(rayPath));
  }  // end of loop over observation samples in the profile (except the topmost obs).

  oops::Log::trace() << name() << "::makeProfileRayPaths - exit "
        << profileRayPaths->numRays() << " ray paths created" << std::endl;
  return profileRayPaths;
}

}  // namespace ufo
