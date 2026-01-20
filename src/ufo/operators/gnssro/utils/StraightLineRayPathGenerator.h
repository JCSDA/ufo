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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_STRAIGHTLINERAYPATHGENERATOR_H_
#define UFO_OPERATORS_GNSSRO_UTILS_STRAIGHTLINERAYPATHGENERATOR_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "ufo/operators/gnssro/utils/GnssroRayPathGenerator.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Subclass of GnssroRayPathGenerator for producing straight-line ray paths.
///
//
class StraightLineRayPathGenerator: public GnssroRayPathGenerator {
 public:
  class NodeCalculator {
   public:
    NodeCalculator(double tpGeomHeight,
                   double localEarthRadius,
                   float azimuth,
                   float tpLat,
                   float tpLon,
                   const GnssroRayPathParameters & params);
    ~NodeCalculator();

    // Accessors
    double tpLat() const {return tpLat_;}
    double tpLon() const {return tpLon_;}
    double azimuth() const {return azimuth_;}
    double edgeDist() const {return edgeDist_;}
    double edgeRadius() const {return edgeRadius_;}
    double edgeGeomHeight() const {return edgeRadius_ - localEarthRadius_;}
    double rightTheta() const {return baseTheta_;}
    double leftTheta() const {return -1.0 * baseTheta_;}

    /// \brief Method to extend ray properties by one more node in
    //  each direction from the tangent point.
    void extendRay();

    /// Returns lon, lat in degrees that is arc-length theta from the
    //  tangent point along the ray path defined by the azimuth angle.
    std::pair<float, float> lonLatAlongRay(double theta);

   private:
    double tpGeomHeight_;
    double localEarthRadius_;
    double radius0_;
    double radius0_sq_;
    double azimuth_;
    double tpLat_;
    double tpLon_;
    GnssroRayPathParameters params_;
    double azimuthRad_;
    double tpLatRad_;
    double tpLonRad_;
    double cosAzimuth_;
    int   signAzimuth_;
    double cosTpLat_;
    double cosTpLon_;
    double sinTpLat_;
    double sinTpLon_;
    double edgeIncrement_;
    double edgeDist_;
    double edgeRadius_;
    double baseTheta_;
  };  // end NodeCalculator nested class

  static constexpr const char GENTYPE[] = "StraightLine";
  static const char * name() {return GENTYPE;}

  //  Constructor/Destructor
  explicit StraightLineRayPathGenerator(
      const ioda::ObsSpace & odb,
      const GnssroRayPathParameters & params);
  ~StraightLineRayPathGenerator() override;

  /// \brief Return a GnssroRayPath object based on observation data read for
  /// the specified GNSSRO profile from the ObsSpace.
  std::unique_ptr<GnssroProfileRayPaths> makeProfileRayPaths(
      const GnssroProfileSlice & profileSlice) override;

 private:
  std::vector<float> obsLat_;      // observed latitude (degrees)
  std::vector<float> obsLon_;      // observed longitude (degrees)
  std::vector<float> obsAlt_;      // observed geometric height (m)
  std::vector<float> obsRoc_;      // observed earth radius of curvature (m)
  std::vector<float> obsUndul_;    // observed geoid undulation (m)
  std::vector<float> obsAzimuth_;  // observed azimuth angle/bearing (deg clockwise from north)
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_STRAIGHTLINERAYPATHGENERATOR_H_
