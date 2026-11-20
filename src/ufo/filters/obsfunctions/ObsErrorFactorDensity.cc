/*
 * (C) Copyright 2026 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorDensity.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

#include "eckit/container/KDTree.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point3.h"
#include "ioda/distribution/Distribution.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorDensity> makerObsErrorFactorDensity_("ObsErrorFactorDensity");

namespace {

/// Position of each observation, and whether it is counted in the density.
struct ObsPositions {
  std::vector<int> counted;     // 1 if the observation is counted, 0 otherwise
  std::vector<int> interval;    // index of the time interval (only set if counted)
  std::vector<float> lat;
  std::vector<float> lon;
  std::vector<float> vcoord;
};

/// KD-tree of observation positions on the unit sphere; the payload is the observation index.
struct TreeTraits {
  using Point = eckit::geometry::Point3;
  using Payload = size_t;
};
using KDTree = eckit::KDTreeMemory<TreeTraits>;

eckit::geometry::Point3 unitSpherePoint(double lat, double lon) {
  const double phi = lat * Constants::deg2rad;
  const double lambda = lon * Constants::deg2rad;
  return eckit::geometry::Point3(std::cos(phi) * std::cos(lambda),
                                 std::cos(phi) * std::sin(lambda),
                                 std::sin(phi));
}

/// Reads the positions of the observations on this task. An observation is counted if it passed
/// QC and has a location, time and vertical coordinate.
ObsPositions readPositions(const ObsFilterData & data,
                           const ObsErrorFactorDensityParameters & options) {
  const size_t nlocs = data.nlocs();
  ObsPositions obs;
  obs.lat.resize(nlocs);
  obs.lon.resize(nlocs);
  obs.vcoord.resize(nlocs);
  std::vector<int> qcflag(nlocs);
  std::vector<util::DateTime> datetime(nlocs);
  data.get(Variable(options.testQCflag.value() + "/" + options.variable.value()), qcflag);
  data.get(Variable("MetaData/latitude"), obs.lat);
  data.get(Variable("MetaData/longitude"), obs.lon);
  data.get(Variable("MetaData/dateTime"), datetime);
  data.get(options.verticalCoordinate.value(), obs.vcoord);

  const float missing = util::missingValue<float>();
  const util::DateTime missingTime = util::missingValue<util::DateTime>();
  const util::DateTime windowStart = data.obsspace().windowStart();
  const int64_t timeMesh = options.timeMesh.value().toSeconds();

  obs.counted.assign(nlocs, 0);
  obs.interval.assign(nlocs, 0);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    obs.counted[jloc] = qcflag[jloc] == 0 && obs.lat[jloc] != missing &&
                        obs.lon[jloc] != missing && obs.vcoord[jloc] != missing &&
                        datetime[jloc] != missingTime;
    if (obs.counted[jloc]) {
      // Like the time window, the intervals (start + k mesh, start + (k + 1) mesh] exclude their
      // beginning and include their end, so an observation at the end of the window is in the last
      // interval rather than in an interval of its own.
      obs.interval[jloc] = ((datetime[jloc] - windowStart).toSeconds() - 1) / timeMesh;
    }
  }
  return obs;
}

/// Gathers the positions from all tasks (each observation appears once).
ObsPositions gatherPositions(const ObsPositions & local, const ioda::Distribution & dist) {
  ObsPositions all(local);
  dist.allGatherv(all.counted);
  dist.allGatherv(all.interval);
  dist.allGatherv(all.lat);
  dist.allGatherv(all.lon);
  dist.allGatherv(all.vcoord);
  return all;
}

/// Builds one KD-tree of the counted observations for each time interval.
std::map<int, KDTree> buildTrees(const ObsPositions & obs) {
  std::map<int, std::vector<KDTree::Value>> points;
  for (size_t jloc = 0; jloc < obs.counted.size(); ++jloc) {
    if (obs.counted[jloc]) {
      points[obs.interval[jloc]].emplace_back(unitSpherePoint(obs.lat[jloc], obs.lon[jloc]), jloc);
    }
  }
  std::map<int, KDTree> trees;
  for (auto & [interval, intervalPoints] : points) {
    trees[interval].build(intervalPoints.begin(), intervalPoints.end());
  }
  return trees;
}

}  // namespace

// -----------------------------------------------------------------------------

ObsErrorFactorDensity::ObsErrorFactorDensity(const eckit::LocalConfiguration & config)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorDensity constructor" << std::endl;
  options_.validateAndDeserialize(config);
  if (options_.timeMesh.value().toSeconds() <= 0) {
    throw eckit::BadValue("ObsErrorFactorDensity: time mesh must be positive", Here());
  }
  invars_ += Variable(options_.testQCflag.value() + "/" + options_.variable.value());
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/longitude");
  invars_ += Variable("MetaData/dateTime");
  invars_ += options_.verticalCoordinate.value();
}

// -----------------------------------------------------------------------------

ObsErrorFactorDensity::~ObsErrorFactorDensity() {
  oops::Log::trace() << "ObsErrorFactorDensity destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorDensity::compute(const ObsFilterData & data,
                                    ioda::ObsDataVector<float> & factor) const {
  oops::Log::trace() << "ObsErrorFactorDensity compute start" << std::endl;

  const ObsPositions local = readPositions(data, options_);
  const ObsPositions all = gatherPositions(local, *data.obsspace().distribution());
  std::map<int, KDTree> trees = buildTrees(all);

  // The counting volume is a cylinder centred on each observation, with the horizontal area of a
  // square of the horizontal mesh size and the height of the vertical mesh size. On the unit
  // sphere, observations within the cylinder's radius are within a chord of 2 sin(angle / 2).
  const double radius = Constants::thousand * options_.horizontalMesh.value() / std::sqrt(M_PI);
  const double angle = std::min(M_PI, radius / Constants::mean_earth_rad_m);
  const double chord = Constants::two * std::sin(Constants::half * angle);
  const double halfHeight = Constants::half * options_.verticalMesh.value();
  const float threshold = options_.densityThreshold.value();

  for (size_t jloc = 0; jloc < local.counted.size(); ++jloc) {
    float density = 0.0f;
    if (local.counted[jloc]) {
      const auto neighbours = trees.at(local.interval[jloc]).findInSphere(
            unitSpherePoint(local.lat[jloc], local.lon[jloc]), chord);
      for (const auto & neighbour : neighbours) {
        if (std::abs(all.vcoord[neighbour.payload()] - local.vcoord[jloc]) <= halfHeight) {
          density += 1.0f;
        }
      }
    }
    const float locFactor = std::sqrt(std::max(1.0f, density / threshold));
    for (size_t ivar = 0; ivar < factor.nvars(); ++ivar) {
      factor[ivar][jloc] = locFactor;
    }
  }
  oops::Log::trace() << "ObsErrorFactorDensity compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorDensity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
