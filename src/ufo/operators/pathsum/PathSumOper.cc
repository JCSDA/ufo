/*
 * (C) Copyright
 *
 * Licensed under the terms of the Apache Licence Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/pathsum/PathSumOper.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"
#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<PathSumOper> makerPathSum_("PathSum");

// -----------------------------------------------------------------------------
PathSumOper::PathSumOper(const ioda::ObsSpace & odb,
                         const Parameters_ & params)
  : ObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb),
    pathType_(params.pathType),
    geovalVar_(oops::Variable(params.geovalVar.value())),
    weightsVar_(params.weightsVar.value() ?
          boost::optional<oops::Variable>(oops::Variable(*params.weightsVar.value()))
          : boost::none),
    weights_(params.weights.value() != boost::none ?
             *params.weights.value() : std::vector<double>()),
    heightRange_(params.heightRange.value() != boost::none ?
                 *params.heightRange.value() : std::vector<double>()),
    interpolateBoundaries_(params.interpolateBoundaries),
    useKmForHeight_(params.useKmForHeight.value()),
    scalingFactor_(params.scalingFactor)
  {
    requiredVars_ += oops::Variables(std::vector<oops::Variable>{geovalVar_});
    if (!weightsVar_ || !heightRange_.empty()) {
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                       oops::Variable("height_wrt_surface")});
    }

    // If weights come from GeoVaLs, require them
    if (weightsVar_) {
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{*weightsVar_});
    }

    oops::Log::trace() << "PathSumOper constructor done." << std::endl;
  }
// -----------------------------------------------------------------------------

PathSumOper::~PathSumOper() {
  oops::Log::trace() << "PathSumOper destructor done" << std::endl;
}

// -----------------------------------------------------------------------------
/// Segement length between two points is returned in [m] by default
/// or [km] if useKmForHeight is true
double PathSumOper::computeSegmentLength(const std::array<double, 3> &p1,
                                         const std::array<double, 3> &p2) const {
  // Apply scaling factor if input heights are in km
  double scale = useKmForHeight_ ? 0.001 : 1.0;                // convert m->km if needed

  // WGS-84 ellipsoid parameters
  const double a_scaled = Constants::semi_major_axis * scale;   // semi-major axis of Earth
                                                                // default [m] or [km] if scaled
  const double b_scaled = Constants::semi_minor_axis * scale;   // semi-minor axis of Earth
                                                                // default [m] or [km] if scaled
  const double ecc2 = Constants::eccentricity_sq;       // eccentricity squared of the ellipsoid

  auto geodeticToECEF = [&](const std::array<double, 3>& p) -> std::array<double, 3> {
    double lat = p[0] * M_PI / 180.0;  // geodetic latitude
    double lon = p[1] * M_PI / 180.0;  // geodetic longitude
    double height   = p[2];            // height wrt Earth surface

    double radius = a_scaled / std::sqrt(1.0 - ecc2 * std::sin(lat)
                                               * std::sin(lat));  // radius of curvature

    double x = (radius + height) * std::cos(lat) * std::cos(lon);
    double y = (radius + height) * std::cos(lat) * std::sin(lon);
    double z = (radius * (1.0 - ecc2) + height) * std::sin(lat);

    return {x, y, z};
  };

  // Convert points to ECEF
  std::array<double, 3> ecef1 = geodeticToECEF(p1);
  std::array<double, 3> ecef2 = geodeticToECEF(p2);

  // Compute Euclidean distance
  double dx = ecef2[0] - ecef1[0];
  double dy = ecef2[1] - ecef1[1];
  double dz = ecef2[2] - ecef1[2];

  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//------------------------------------------------------------------------------
/// Check if a value is valid (not missing, NaN, or extreme)
inline bool isValidValue(double v) {
  const double missingValueDouble = util::missingValue<double>();
  const double missingValueFloat = util::missingValue<float>();

  // Reject NaN or Inf
  if (!std::isfinite(v)) return false;

  // Reject explicit UFO missing values
  if (v == static_cast<double>(missingValueFloat) ||
      v == missingValueDouble)
    return false;

  return true;
}

//------------------------------------------------------------------------------
/// TrapezoidalWeight = 0.5 * (segment ending at lev) + 0.5 * (segment starting at lev)
double PathSumOper::computeTrapezoidalWeight(std::size_t lev,
                                             const std::vector<double> &heights,
                                             float lat, float lon,
                                             std::size_t nlevs) const {
  const double missing = util::missingValue<double>();

  if (heights.empty() || lev >= nlevs) {
    return missing;
  }

  double prevSegmentLength = missing;
  double nextSegmentLength = missing;

  // Segment ending at lev (from lev-1 to lev)
  if (lev > 0 && isValidValue(heights[lev-1]) && isValidValue(heights[lev])) {
    std::array<double, 3> p0 = {lat, lon, heights[lev-1]};
    std::array<double, 3> p1 = {lat, lon, heights[lev]};
    prevSegmentLength = computeSegmentLength(p0, p1);
  }

  // Segment starting at lev (from lev to lev+1)
  if (lev < nlevs - 1 && isValidValue(heights[lev]) && isValidValue(heights[lev+1])) {
    std::array<double, 3> p0 = {lat, lon, heights[lev]};
    std::array<double, 3> p1 = {lat, lon, heights[lev+1]};
    nextSegmentLength = computeSegmentLength(p0, p1);
  }

  // Trapezoidal weight: half of each adjacent segment
  if (isValidValue(prevSegmentLength) && isValidValue(nextSegmentLength)) {
    return 0.5 * prevSegmentLength + 0.5 * nextSegmentLength;
  } else if (isValidValue(nextSegmentLength)) {
    return 0.5 * nextSegmentLength;  // First level: only next segment
  } else if (isValidValue(prevSegmentLength)) {
    return 0.5 * prevSegmentLength;  // Last level: only prev segment
  }

  return missing;
}

//------------------------------------------------------------------------------
/// Linear interpolation to get value at point x between two points x0 (with
/// value y0) and x1 (with value y1)
inline double linearInterpolate(double x, double x0, double x1, double y0, double y1) {
  if (x1 == x0) return y0;  // avoid div-by-zero
  double slope = (y1 - y0) / (x1 - x0);
  double yinterp = y0 + (x - x0) * slope;
  return yinterp;
}

// -----------------------------------------------------------------------------
void PathSumOper::simulateObs(const GeoVaLs & geovals,
                              ioda::ObsVector & hofx,
                              ObsDiagnostics &, const QCFlags_t &) const {
// PathSum: Integration of a geoval variable along a vertical or any slant path
// defined by latitude, longitude and height

  oops::Log::trace() << "PathSumOper::simulateObs start" << std::endl;

  const std::size_t nlocs = geovals.nlocs();
  const std::size_t nlevs = geovals.nlevs(geovalVar_);
  ASSERT(nlocs == hofx.nlocs());
  ASSERT(nlocs == odb_.nlocs());

  const double missing = util::missingValue<double>();

  // Lat/lon from ODB only needed if computing weights(wts)
  std::vector<float> lats, lons;
  if (!weightsVar_) {
    odb_.get_db("MetaData", "latitude", lats);
    odb_.get_db("MetaData", "longitude", lons);
    ASSERT(lats.size() == nlocs && lons.size() == nlocs);
  }

  if (pathType_ == "vertical") {
    // Case 1: Vertical integration

    // Loop over locations
    for (std::size_t loc = 0; loc < nlocs; ++loc) {
      hofx[loc] = missing;

      // Extract profile for main variable
      std::vector<double> valueProfile(nlevs);
      geovals.getAtLocation(valueProfile, geovalVar_, loc);

      // Extract weights for integration computation if WeightsVar is defined in yaml
      std::vector<double> wts(nlevs, missing);  // default = missing, dimension nlevs
      if (weightsVar_) {
        geovals.getAtLocation(wts, *weightsVar_, loc);
      }

      // Apply height range restriction and interpolate boundaries if needed
      const bool useHeightRange = !heightRange_.empty();
      const double hmin = useHeightRange ? heightRange_[0] : std::numeric_limits<double>::lowest();
      const double hmax = useHeightRange ? heightRange_[1] : std::numeric_limits<double>::max();

      // Get vertical coordinate if no weights are defined or height range is defined in yaml
      std::vector<double> heights;
      if ((!weightsVar_ && weights_.empty()) || useHeightRange) {
        std::vector<double> heightProfile(nlevs);
        geovals.getAtLocation(heightProfile, oops::Variable{"height_wrt_surface"}, loc);
        if (heightProfile.empty()) {
          oops::Log::warning() << "height_wrt_surface is empty," <<
                                  " therefore height range can not be applied" << std::endl;
          continue;  // safety check
        }
        heights = heightProfile;  // Copy to heights for use in interpolation and integration
      }

      // Create working copies so original GeoVaLs remain untouched
      std::vector<double> vals = valueProfile;

      if (useHeightRange && interpolateBoundaries_ && heights.size() > 1) {
        // Ensure increasing order. Note: heights should be monotonic
        if (heights.front() > heights.back()) {
          std::reverse(heights.begin(), heights.end());
          std::reverse(vals.begin(), vals.end());
        }

        // Interpolate upper/lower boundary if needed
        bool lowerBoundaryInterpolated = false;
        bool upperBoundaryInterpolated = false;
        bool lowerBoundaryNeeded = false;
        bool upperBoundaryNeeded = false;
        for (std::size_t lev = 1; lev < nlevs; ++lev) {
          if (isValidValue(heights[lev-1]) && heights[lev-1] < hmin) {
            lowerBoundaryNeeded = true;
            if (isValidValue(heights[lev]) && heights[lev] >= hmin) {
              double interpVal = linearInterpolate(hmin,
                                                   heights[lev-1], heights[lev],
                                                   vals[lev-1], vals[lev]);
              heights[lev-1] = hmin;
              vals[lev-1] = interpVal;
              lowerBoundaryInterpolated = true;
              break;
            }
          }
          if (isValidValue(heights[lev]) && heights[lev] >= hmax) {
            upperBoundaryNeeded = true;
            if (isValidValue(heights[lev-1]) && heights[lev-1] < hmax) {
              double interpVal = linearInterpolate(hmax,
                                                   heights[lev-1], heights[lev],
                                                   vals[lev-1], vals[lev]);
              heights[lev] = hmax;
              vals[lev] = interpVal;
              upperBoundaryInterpolated = true;
              break;
            }
          }
        }
        if (lowerBoundaryNeeded && !lowerBoundaryInterpolated) {
            oops::Log::warning() << "Could not interpolate lower boundary at hmin=" << hmin
                                 << " for location " << loc
                                 << ": no valid level pair found spanning hmin" << std::endl;
        }
        if (upperBoundaryNeeded && !upperBoundaryInterpolated) {
           oops::Log::warning() << "Could not interpolate upper boundary at hmax=" << hmax
                                 << " for location " << loc
                                 << ": no valid level pair found spanning hmax" << std::endl;
        }
      }

      // Prepare for vertical path integration
      // Loop over levels (nlevs dimension, same as vals)
      double sum = 0.0;
      bool anyValid = false;

      for (std::size_t lev = 0; lev < nlevs; ++lev) {
        // Check if this level is within height range (if defined)
        bool withinHeightRange = true;
        if (useHeightRange) {
          if (!heights.empty() && isValidValue(heights[lev])) {
            withinHeightRange = (heights[lev] >= hmin && heights[lev] <= hmax);
          } else {
            withinHeightRange = false;
          }
        }

        if (!withinHeightRange) {
          oops::Log::trace() << "Beyond height range, skip level " << heights[lev] << std::endl;
          continue;
        }

        // ---- Determine weight (wt) depending on configuration ----
        double wt = missing;
        if (weightsVar_) {
          // option 1: read from GeoVaLs (dimension nlevs)
          if (lev < wts.size()) {
            wt = wts[lev];
            oops::Log::trace() << "wts from geoval" << std::endl;
          }
        } else if (!weights_.empty()) {
          if (lev < weights_.size()) {
            wt = weights_[lev];
            oops::Log::trace() << "wts from yaml" << std::endl;
          }
        } else {
          // option 3: compute trapezoidal weight based on segment length
          wt = computeTrapezoidalWeight(lev, heights, lats[loc], lons[loc], nlevs);
          oops::Log::trace() << "wts from operator" << std::endl;
        }

        // Compute sum using weight * value
        if (isValidValue(vals[lev]) && isValidValue(wt)) {
          sum += wt * vals[lev];
          anyValid = true;
        }
      }

      // Store result in hofx
      if (anyValid) {
        // convert to certain unit if needed. scalingFactor default is 1.0
        hofx[loc] = sum * scalingFactor_;
      }
    }
  } else if (pathType_ == "slant") {
      oops::Log::warning() << "[PathSumOper] Slant path computation is currently a placeholder "
                           << "and will be implemented in a future update. "
                           << "Exiting now to prevent use of incomplete results." << std::endl;

      throw eckit::NotImplemented
                ("Slant path computation not yet implemented in PathSumOper.", Here());
      return;
  }  // End of vertical or slant path

  oops::Log::trace() << "PathSumOper::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

  void PathSumOper::print(std::ostream & os) const {
  os << "PathSumOper::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace ufo
