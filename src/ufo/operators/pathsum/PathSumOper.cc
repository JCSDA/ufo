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
double PathSumOper::computeSegmentLength(const std::array<double, 3> &p1,
                                         const std::array<double, 3> &p2) const {
// Distance between two points is retruned is in [m] by default
// or [km] if useKmForHeight is true

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
  std::vector<double> lats, lons;
  if (!weightsVar_) {
    odb_.get_db("MetaData", "latitude", lats);
    odb_.get_db("MetaData", "longitude", lons);
    ASSERT(lats.size() == nlocs && lons.size() == nlocs);
  }

  if (pathType_ == "vertical") {
    // Case 1: Vertical integration
    // Data structure requirement:
    // for each location defined by latitude(Location), longitude(Location), there is
    // a vertical model profile defined by height_wrt_surface(Location, nlevs)

    // Loop over locations
    for (std::size_t loc = 0; loc < nlocs; ++loc) {
      hofx[loc] = missing;

      // Extract profile for main variable
      std::vector<double> vals(nlevs);
      geovals.getAtLocation(vals, geovalVar_, loc);

      // Get vertical coordinate if no weights are defined or height range is defined in yaml
      std::vector<double> heights(nlevs);
      if ((!weightsVar_ && weights_.empty()) || !heightRange_.empty()) {
        geovals.getAtLocation(heights, oops::Variable{"height_wrt_surface"}, loc);
        if (heights.empty()) {
          oops::Log::warning() << "height_wrt_surface is empty," <<
                                  " therefore height range can not be applied" << std::endl;
          continue;  // safety check
        }
      }

      // Extract weights for integration computation if WeightsVar is defined in yaml
      std::vector<double> wts(nlevs-1, missing);  // default = missing
      if (weightsVar_) {
        geovals.getAtLocation(wts, *weightsVar_, loc);
      }

      // Prepare for vertical path integration
      double sum = 0.0;
      bool anyValid = false;

      const bool useHeightRange = !heightRange_.empty();
      const double hmin = useHeightRange ? heightRange_[0] : std::numeric_limits<double>::lowest();
      const double hmax = useHeightRange ? heightRange_[1] : std::numeric_limits<double>::max();

      for (std::size_t lev = 1; lev < nlevs; ++lev) {
        double wt = missing;

        // Determine whether this level pair should be considered based on height range (if defined)
        bool withinHeightRange = true;
        if (!heightRange_.empty()) {
          // If height is available and height range is define in yaml,
          // only consider levels within the height range
          withinHeightRange =
            (heights[lev]   >= hmin && heights[lev]   < hmax &&
             heights[lev-1] >= hmin && heights[lev-1] < hmax);
        }

        if (!withinHeightRange) {
          oops::Log::trace() << "Beyond height range, skip level " << lev << std::endl;
          continue;  // skip levels outside the height range
        }

        // ---- Determine weight (wt) depending on configuration ----
        if (weightsVar_) {
          // option 1: read from GeoVaLs
          wt = wts[lev-1];  // pre-fetched earlier
        } else if (!weights_.empty()) {
          // option 2: yaml-defined static weights
          if (lev - 1 < weights_.size())
            wt = weights_[lev - 1];
        } else {
          // option 3: compute geometrically (requires height + lat/lon)
          if (!heights.empty() &&
            isValidValue(heights[lev-1]) && isValidValue(heights[lev])) {
            std::array<double, 3> p0 = {lats[loc], lons[loc], heights[lev-1]};
            std::array<double, 3> p1 = {lats[loc], lons[loc], heights[lev]};
            wt = computeSegmentLength(p0, p1);  // m or km depending on flag
          }
        }

        // Compute trapezoidal sum if all valid
        if (isValidValue(vals[lev]) && isValidValue(vals[lev-1]) && isValidValue(wt)) {
          sum += 0.5f * (vals[lev] + vals[lev-1]) * wt;
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
