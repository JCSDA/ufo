/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/Sphere.h"
#include "ioda/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/sqr.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

namespace {
/// \brief Returns the square of the distance between the two \p Point arguments.
float distance2(const TrackCheckUtils::Point &a, const TrackCheckUtils::Point &b) {
  float sum = 0;
  for (size_t i = 0; i < a.size(); ++i)
    sum += util::sqr(a[i] - b[i]);
  return sum;
}

TrackCheckUtils::Point pointFromLatLon(float latitude, float longitude) {
  // This local copy is needed because convertSphericalToCartesian takes the first parameter by
  // reference, but Constants::mean_earth_rad has no out-of-line definition.
  const double meanEarthRadius = Constants::mean_earth_rad;
  eckit::geometry::Point3 eckitPoint;
  eckit::geometry::Sphere::convertSphericalToCartesian(
        meanEarthRadius, eckit::geometry::Point2(longitude, latitude), eckitPoint);
  TrackCheckUtils::Point Point;
  std::copy(eckitPoint.begin(), eckitPoint.end(), Point.begin());
  return Point;
}
}  // namespace

/// \p Returns the distance between the two cartesian-mapped \p Point arguments
float TrackCheckUtils::distance(const Point &a, const Point &b) {
  return std::sqrt(distance2(a, b));
}

ObsAccessor TrackCheckUtils::createObsAccessor(const boost::optional<Variable> &stationIdVariable,
                                               const ioda::ObsSpace &obsdb) {
  if (stationIdVariable != boost::none) {
    return ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
          obsdb, *stationIdVariable);
  } else if (!obsdb.obs_group_var().empty()) {
    // Assume observations were grouped into records by station IDs
    return ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(obsdb);
  } else {
    // Observations were not grouped into records and no station ID variable was provided.
    // Assume all observations were taken by the same station.
    return ObsAccessor::toAllObservations(obsdb);
  }
}

void TrackCheckUtils::sortTracksChronologically(const std::vector<size_t> &validObsIds,
                                                const ObsAccessor &obsAccessor,
                                                RecursiveSplitter &splitter) {
  const std::vector<util::DateTime> times = obsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "datetime");
  splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndexA, size_t obsIndexB)
  { return times[validObsIds[obsIndexA]] < times[validObsIds[obsIndexB]]; });
}

TrackCheckUtils::ObsGroupLocationTimes
  TrackCheckUtils::collectObservationsLocations(const ObsAccessor &obsAccessor) {
  ObsGroupLocationTimes locationTimes;

  locationTimes.latitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
  locationTimes.longitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "longitude");
  locationTimes.datetimes = obsAccessor.getDateTimeVariableFromObsSpace("MetaData", "datetime");

  return locationTimes;
}

TrackCheckUtils::ObsLocationTime::ObsLocationTime(float latitude, float longitude,
                                                                  const util::DateTime &time)
  :  location_(pointFromLatLon(latitude, longitude)), time_(time)
{}

float TrackCheckUtils::CheckCounter::failedChecksFraction() const {
  return numChecks_ != 0 ? static_cast<float>(numFailedChecks_) / numChecks_ : 0.0f;
}

TrackCheckUtils::CheckCounter::CheckCounter()
  :  numFailedChecks_(0), numChecks_(0)
{}

void TrackCheckUtils::CheckCounter::registerCheckResult(const CheckResult &result) {
  if (result != CheckResult::SKIPPED)
    ++numChecks_;
  if (result == CheckResult::FAILED)
    ++numFailedChecks_;
  assert(numChecks_ >= 0);
  assert(numFailedChecks_ >= 0);
  assert(numFailedChecks_ <= numChecks_);
}

void TrackCheckUtils::CheckCounter::unregisterCheckResult(const CheckResult &result) {
  if (result != CheckResult::SKIPPED)
    --numChecks_;
  if (result == CheckResult::FAILED)
    --numFailedChecks_;
  assert(numChecks_ >= 0);
  assert(numFailedChecks_ >= 0);
  assert(numFailedChecks_ <= numChecks_);
}

}  // namespace ufo
