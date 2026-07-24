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
    scalingFactor_(params.scalingFactor),
    pathPointLatVar_(oops::Variable(params.pathPointLatVar.value())),
    pathPointLonVar_(oops::Variable(params.pathPointLonVar.value())),
    pathPointHeightVar_(oops::Variable(params.pathPointHeightVar.value()))
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

    if (pathType_ == PathType::SLANT && !weightsVar_ && weights_.empty()) {
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{*pathPointLatVar_});
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{*pathPointLonVar_});
    }
    if ((!weightsVar_ && weights_.empty()) || !heightRange_.empty() ||
         interpolateBoundaries_) {
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{*pathPointHeightVar_});
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
/// NOTE: Assume the latitudes from geovals are geodetic latitudes
double PathSumOper::computeSegmentLength(const std::array<double, 3> &p1,
                                         const std::array<double, 3> &p2) const {
  // Apply scaling factor if input heights are in km
  double scale = useKmForHeight_ ? 0.001 : 1.0;

  // WGS-84 ellipsoid parameters
  const double a_scaled = Constants::semi_major_axis * scale;   // semi-major axis of Earth
                                                                // default [m] or [km] if scaled
  const double ecc2 = Constants::eccentricity_sq;       // eccentricity squared of the ellipsoid

  auto geodeticToECEF = [&](const std::array<double, 3>& p) -> std::array<double, 3> {
    double lat = p[0] * M_PI / 180.0;  // geodetic latitude
    double lon = p[1] * M_PI / 180.0;  // geodetic longitude
    double height   = p[2];            // geodetic height

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
double PathSumOper::computeTrapezoidalWeight(std::size_t ipoint,
                                             const PathGeometry &geom) const {
  const double missing = util::missingValue<double>();
  const std::size_t npoint = geom.height.size();

  if (npoint == 0 || ipoint >= npoint) return missing;

  double prevSegmentLength = missing;
  double nextSegmentLength = missing;

  // Segment ending at ipoint (from ipoint-1 to ipoint)
  if (ipoint > 0 && isValidValue(geom.height[ipoint-1]) &&
                    isValidValue(geom.height[ipoint])) {
    std::array<double, 3> p0 = {geom.lat[ipoint-1], geom.lon[ipoint-1], geom.height[ipoint-1]};
    std::array<double, 3> p1 = {geom.lat[ipoint], geom.lon[ipoint], geom.height[ipoint]};
    prevSegmentLength = computeSegmentLength(p0, p1);
  }

  // Segment starting at ipoint (from ipoint to ipoint+1)
  if (ipoint < npoint - 1 && isValidValue(geom.height[ipoint]) &&
                              isValidValue(geom.height[ipoint+1])) {
    std::array<double, 3> p0 = {geom.lat[ipoint], geom.lon[ipoint], geom.height[ipoint]};
    std::array<double, 3> p1 = {geom.lat[ipoint+1], geom.lon[ipoint+1], geom.height[ipoint+1]};
    nextSegmentLength = computeSegmentLength(p0, p1);
  }

  // Trapezoidal weight: half of each adjacent segment
  if (isValidValue(prevSegmentLength) && isValidValue(nextSegmentLength)) {
    return 0.5 * prevSegmentLength + 0.5 * nextSegmentLength;
  } else if (isValidValue(nextSegmentLength)) {
    return 0.5 * nextSegmentLength;  // First ipoint: only next segment
  } else if (isValidValue(prevSegmentLength)) {
    return 0.5 * prevSegmentLength;  // Last ipoint: only prev segment
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
    oops::Log::warning() << "Warning:pathSumOper::interpolateBoundaries: " <<
          "interpolateBoundaries called with less than 2 heights" << std::endl;
    return;
  }

  const std::size_t nlevs = heights.size();
  HeightOrder order = detectHeightOrder(heights);
  if (order == HeightOrder::UNKNOWN) {
    oops::Log::warning()
          << "Warning:pathSumOper::interpolateBoundaries: "
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
    oops::Log::warning() << "pathSumOper::interpolateBoundaries: "
                       << "Could not interpolate lower boundary at hmin=" << hmin << std::endl;
  if (upperNeeded && !upperInterpolated)
    oops::Log::warning() << "pathSumOper::interpolateBoundaries: "
                       << "Could not interpolate upper boundary at hmax=" << hmax << std::endl;

  // ---- mask values outside height range ----
  const double missing = util::missingValue<double>();

  for (std::size_t lev = 0; lev < heights.size(); ++lev) {
    if (!isValidValue(heights[lev]) ||
        heights[lev] < hmin ||
        heights[lev] > hmax) {
      vals[lev] = missing;
    }
  }
}

// -----------------------------------------------------------------------------
PathSumOper::PathGeometry
PathSumOper::buildGeometry(std::size_t loc,
                           const GeoVaLs & geovals,
                           const std::vector<float> &lats,
                           const std::vector<float> &lons,
                           bool needHeights,
                           bool needPathLatLon) const {
  PathGeometry geom;

  const std::size_t npoints = geovals.nlevs(geovalVar_);
  const float missingValueFloat = util::missingValue<float>();
  const double missingValueDouble = util::missingValue<double>();

  geom.lat.assign(npoints, missingValueFloat);
  geom.lon.assign(npoints, missingValueFloat);
  geom.height.assign(npoints, missingValueDouble);

  if (needHeights) {
    try {
      geovals.getAtLocation(geom.height, oops::Variable{*pathPointHeightVar_}, loc);
    } catch (const std::exception &error) {
      oops::Log::warning()
        << "Warning:PathSumOper: failed to read " << *pathPointHeightVar_
        << " at loc=" << loc << " : " << error.what() << std::endl;
      throw;
    }
  }

  if (pathType_ == PathType::VERTICAL) {
    geom.lat.assign(npoints, lats[loc]);
    geom.lon.assign(npoints, lons[loc]);
  } else {
    if (needPathLatLon) {
      geovals.getAtLocation(geom.lat, oops::Variable{*pathPointLatVar_}, loc);
      geovals.getAtLocation(geom.lon, oops::Variable{*pathPointLonVar_}, loc);
    }
  }

  return geom;
}

// -----------------------------------------------------------------------------
void PathSumOper::buildWeights(std::vector<double> &wts,
                               const PathGeometry &geom,
                               std::size_t loc,
                               const GeoVaLs &geovals,
                               bool useHeightRange,
                               const double hmin, const double hmax,
                               const std::size_t npoint) const {
  const double missing = util::missingValue<double>();
  std::vector<double> wtsValue(npoint, missing);

  if (weightsVar_) {
    geovals.getAtLocation(wtsValue, *weightsVar_, loc);
    oops::Log::trace() << "PathSumOper: using weights from GeoVaLs" << std::endl;
  } else if (!weights_.empty()) {
    std::size_t ncopy = std::min(weights_.size(), npoint);
    std::copy(weights_.begin(), weights_.begin()+ncopy, wtsValue.begin());
    oops::Log::trace() << "PathSumOper: using weights from YAML" << std::endl;
  } else {
    oops::Log::trace() << "PathSumOper: using operator-computed weights" << std::endl;
  }

  for (std::size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    bool withinHeightRange = !useHeightRange ||
                             (!geom.height.empty() && isValidValue(geom.height[ipoint]) &&
                             geom.height[ipoint] >= hmin && geom.height[ipoint] <= hmax);

    if (!withinHeightRange) continue;

    if (weightsVar_ || !weights_.empty()) {
      wts[ipoint]= wtsValue[ipoint];
    } else {
      wts[ipoint] = computeTrapezoidalWeight(ipoint, geom);
    }
  }
}

// -----------------------------------------------------------------------------
double PathSumOper::integrate(const std::vector<double> &vals,
                              const std::vector<double> &wts) const {
  const double missing = util::missingValue<double>();

  double sum = 0.0;
  bool anyValid = false;

  for (std::size_t ipoint = 0; ipoint < vals.size(); ++ipoint) {
    if (!isValidValue(vals[ipoint]) || !isValidValue(wts[ipoint])) continue;
    sum += vals[ipoint] * wts[ipoint];
    anyValid = true;
  }

  return anyValid ? sum : missing;
}
// -----------------------------------------------------------------------------
void PathSumOper::simulateObs(const GeoVaLs & geovals,
                              ioda::ObsVector & hofx,
                              ObsDiagnostics &, const QCFlags_t &) const {
  const std::size_t nlocs = geovals.nlocs();  // corresponds to observation locations
  const std::size_t npoints = geovals.nlevs(geovalVar_);   // vertical levels or slant path points
  const double missing = util::missingValue<double>();
  const bool useHeightRange = !heightRange_.empty();

  std::vector<float> lats, lons;
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);

  ASSERT(lats.size() == nlocs);
  ASSERT(lons.size() == nlocs);

  // Decide if heights are required
  const bool needHeights = (!weightsVar_ && weights_.empty()) ||
                           useHeightRange || interpolateBoundaries_;

  // Decide if path point Lat/Lon are required
  const bool needPathLatLon = (pathType_ == PathType::SLANT && !weightsVar_ && weights_.empty());

  // Height range is defined via YAML or infinite
  const double hmin = useHeightRange ? heightRange_[0] : std::numeric_limits<double>::lowest();
  const double hmax = useHeightRange ? heightRange_[1] : std::numeric_limits<double>::max();

  for (std::size_t loc = 0; loc < nlocs; ++loc) {
    hofx[loc] = missing;

    // Build unified geometry for both vertical and slant paths
    PathGeometry geom = buildGeometry(loc, geovals, lats, lons, needHeights, needPathLatLon);

    // Get geoval values
    std::vector<double> vals(npoints);
    geovals.getAtLocation(vals, geovalVar_, loc);

    if (useHeightRange && interpolateBoundaries_) {
      if (pathType_ == PathType::VERTICAL) {
        interpolateBoundaries(geom.height, vals, hmin, hmax);
      } else {
        oops::Log::warning() << "Warning:PathSumOper:: interpolateBoundaries only implemented "
                           << "for path type: vertical. Ignored the yaml option." << std::endl;
      }
    }

    // Build weights
    std::vector<double> wts(npoints, missing);
    buildWeights(wts, geom, loc, geovals, useHeightRange, hmin, hmax, npoints);

    // Integrate
    double integral = integrate(vals, wts);

    hofx[loc] = isValidValue(integral)
                ? integral * scalingFactor_
                : missing;
  }
  oops::Log::trace() << "PathSumOper::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

  void PathSumOper::print(std::ostream & os) const {
  os << "PathSumOper::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace ufo
