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
    std::string ptype = params.pathType.value();
    std::transform(ptype.begin(), ptype.end(), ptype.begin(), ::tolower);
    if (ptype == "vertical") {
      pathType_ = PathType::VERTICAL;
    } else if (ptype == "slant") {
      pathType_ = PathType::SLANT;
    } else {
      throw eckit::BadParameter(
        "Unknown path type: " + ptype, Here());
    }
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
/// Assumes heights are ordered monotonicly, either increasing or decreasing
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
// Detect height order
enum class HeightOrder {INCREASING, DECREASING, UNKNOWN};
HeightOrder detectHeightOrder(const std::vector<double> &heights) {
  double prev = util::missingValue<double>();
  for (double hcurr : heights) {
    if (!isValidValue(hcurr)) continue;

    if (isValidValue(prev)) {
      if (hcurr > prev) return HeightOrder::INCREASING;
      if (hcurr < prev) return HeightOrder::DECREASING;
      // if equal, keep scanning
    }
    prev = hcurr;
  }
  return HeightOrder::UNKNOWN;
}

// -----------------------------------------------------------------------------
// Interpolation function (order-agnostic, missing-value safe)
void PathSumOper::interpolateBoundaries(std::vector<double> &heights,
                                        std::vector<double> &vals,
                                        double hmin, double hmax) const
{
  if (heights.size() < 2) {
    oops::Log::warning() << "pathSumOper::interpolateBoundaries: " <<
          "interpolateBoundaries called with less than 2 heights" << std::endl;
    return;
  }

  const std::size_t nlevs = heights.size();
  HeightOrder order = detectHeightOrder(heights);
  if (order == HeightOrder::UNKNOWN) {
    oops::Log::warning()
          << "pathSumOper::interpolateBoundaries: "
          << "Height order could not be determined (missing or non-monotonic heights). "
          << "Boundary interpolation will be skipped for such profiles. "
          << std::endl;
    return;
  }

  bool lowerInterpolated = false, upperInterpolated = false;
  bool lowerNeeded = false, upperNeeded = false;

  for (std::size_t lev = 1; lev < nlevs; ++lev) {  // Increasing
    if (order == HeightOrder::INCREASING) {
      // Lower boundary
      if (isValidValue(heights[lev-1]) && heights[lev-1] < hmin) {
        lowerNeeded = true;
        if (isValidValue(heights[lev]) && heights[lev] >= hmin) {
          vals[lev-1] = linearInterpolate(hmin,
                        heights[lev-1], heights[lev], vals[lev-1], vals[lev]);
          heights[lev-1] = hmin;
          lowerInterpolated = true;
          break;
        }
      }
      // Upper boundary
      if (isValidValue(heights[lev]) && heights[lev] >= hmax) {
        upperNeeded = true;
        if (isValidValue(heights[lev-1]) && heights[lev-1] < hmax) {
          vals[lev] = linearInterpolate(hmax,
                      heights[lev-1], heights[lev], vals[lev-1], vals[lev]);
          heights[lev] = hmax;
          upperInterpolated = true;
          break;
        }
      }
    } else {  // Decreasing
      // Lower boundary
      if (isValidValue(heights[lev-1]) && heights[lev-1] > hmin) {
        lowerNeeded = true;
        if (isValidValue(heights[lev]) && heights[lev] <= hmin) {
          vals[lev-1] = linearInterpolate(hmin,
                        heights[lev-1], heights[lev], vals[lev-1], vals[lev]);
          heights[lev-1] = hmin;
          lowerInterpolated = true;
          break;
        }
      }
      // Upper boundary
      if (isValidValue(heights[lev]) && heights[lev] <= hmax) {
        upperNeeded = true;
        if (isValidValue(heights[lev-1]) && heights[lev-1] > hmax) {
          vals[lev] = linearInterpolate(hmax,
                      heights[lev-1], heights[lev], vals[lev-1], vals[lev]);
          heights[lev] = hmax;
          upperInterpolated = true;
          break;
        }
      }
    }
  }

  if (lowerNeeded && !lowerInterpolated)
    oops::Log::debug() << "pathSumOper::interpolateBoundaries: "
                       << "Could not interpolate lower boundary at hmin=" << hmin << std::endl;
  if (upperNeeded && !upperInterpolated)
    oops::Log::debug() << "pathSumOper::interpolateBoundaries: "
                       << "Could not interpolate upper boundary at hmax=" << hmax << std::endl;
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

  if (weightsVar_) {
    oops::Log::debug() << "PathSumOper: using weights from GeoVaLs" << std::endl;
  } else if (!weights_.empty()) {
    oops::Log::debug() << "PathSumOper: using weights from YAML" << std::endl;
  } else {
    oops::Log::debug() << "PathSumOper: using operator-computed weights" << std::endl;
  }

  // Apply height range restriction and interpolate boundaries if needed
  const bool useHeightRange = !heightRange_.empty();
  const double hmin = useHeightRange ? heightRange_[0] : std::numeric_limits<double>::lowest();
  const double hmax = useHeightRange ? heightRange_[1] : std::numeric_limits<double>::max();

  // Vertical integration
  if (pathType_ == PathType::VERTICAL) {
    // Loop over locations
    for (std::size_t loc = 0; loc < nlocs; ++loc) {
      hofx[loc] = missing;

      // Extract profile for main variable
      std::vector<double> vals(nlevs);
      geovals.getAtLocation(vals, geovalVar_, loc);

      // Extract weights for integration computation if WeightsVar is defined in yaml
      std::vector<double> wts(nlevs, missing);  // default = missing, dimension nlevs
      // weights from geovals
      if (weightsVar_) {
        geovals.getAtLocation(wts, *weightsVar_, loc);
      // weights from yaml
      } else if (!weights_.empty()) {
        const std::size_t ncopy_size = std::min(weights_.size(), nlevs);
        std::copy(weights_.begin(), weights_.begin() + ncopy_size, wts.begin());
      }

      // Get vertical coordinate if no weights are defined or height range is defined in yaml
      // or interpolateBoundaries is true
      const bool needHeights = useHeightRange || interpolateBoundaries_ ||
                               (!weightsVar_ && weights_.empty());
      std::vector<double> heights(nlevs);
      if (needHeights) {
        geovals.getAtLocation(heights, oops::Variable{"height_wrt_surface"}, loc);

        if (heights.empty()) {
          oops::Log::debug() << "pathSumOper: height_wrt_surface is empty, PathSum skipped at loc"
                            << loc << std::endl;
            continue;  // safety check
        }
      }

      // Interpolate boundaries
      if (useHeightRange && interpolateBoundaries_ && heights.size() > 1) {
          interpolateBoundaries(heights, vals, hmin, hmax);
      }

      // Prepare for vertical path integration
      // Loop over levels (nlevs dimension, same as vals)
      double sum = 0.0;
      bool anyValid = false;

      for (std::size_t lev = 0; lev < nlevs; ++lev) {
        // Check if this level is within height range (if defined)
        bool withinHeightRange = !useHeightRange ||
                                (!heights.empty() && isValidValue(heights[lev]) &&
                                 heights[lev] >= hmin && heights[lev] <= hmax);

        if (!withinHeightRange || !isValidValue(vals[lev])) continue;

        // ---- Determine weight (wt) depending on configuration ----
        double wt = missing;
        if (lev < wts.size() && isValidValue(wts[lev])) {
          wt = wts[lev];
        } else {
          wt = computeTrapezoidalWeight(lev, heights, lats[loc], lons[loc], nlevs);
        }

        // Compute sum using weight * value
        if (!isValidValue(wt)) continue;

        sum += wt * vals[lev];
        anyValid = true;
      }

      // Store result in hofx and convert to certain unit if needed. scalingFactor default is 1.0
      if (anyValid) hofx[loc] = sum * scalingFactor_;
    }
  }

  // Slant path integration
  if (pathType_ == PathType::SLANT) {
      oops::Log::warning() << "PathSumOper: "
                           << "Slant path computation is currently a placeholder "
                           << "and will be implemented in a future update. "
                           << "Exiting now to prevent use of incomplete results." << std::endl;

      throw eckit::NotImplemented
                ("Slant path computation not yet implemented in PathSumOper.", Here());
  }  // End of vertical or slant path

  oops::Log::trace() << "PathSumOper::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

  void PathSumOper::print(std::ostream & os) const {
  os << "PathSumOper::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace ufo
