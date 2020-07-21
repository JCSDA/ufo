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
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/sqr.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"
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

/// Return the vector of elements of \p categories with indices \p validObsIds.
template <typename T>
std::vector<T> getValidObservationCategories(const std::vector<T> &categories,
                                             const std::vector<size_t> validObsIds)  {
  std::vector<T> validObsCategories(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    validObsCategories[validObsIndex] = categories[validObsIds[validObsIndex]];
  }
  return validObsCategories;
}

std::vector<size_t> TrackCheckUtils::getValidObservationIds(
    const std::vector<bool> &apply, const boost::shared_ptr<ioda::ObsDataVector<int>> &flags) {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && (*(flags))[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsId);
  return validObsIds;
}

void TrackCheckUtils::groupObservationsByStation(const std::vector<size_t> &validObsIds,
                                                 RecursiveSplitter &splitter,
                                                 const eckit::Configuration &config,
                                                 const ioda::ObsSpace &obsdb) {
  std::unique_ptr<TrackCheckUtilsParameters> baseOptions_;
  baseOptions_.reset(new TrackCheckUtilsParameters());
  baseOptions_->deserialize(config);
  if (baseOptions_->stationIdVariable.value() == boost::none) {
    if (obsdb.obs_group_var().empty()) {
      // Observations were not grouped into records.
      // Assume all observations were taken during the same station.
      return;
    } else {
      groupObservationsByRecordNumber(validObsIds, splitter, obsdb);
    }
  } else {
    groupObservationsByVariable(*baseOptions_->stationIdVariable.value(),
                                validObsIds, splitter, obsdb);
  }
}
void TrackCheckUtils::groupObservationsByRecordNumber(const std::vector<size_t> &validObsIds,
                                                      RecursiveSplitter &splitter,
                                                      const ioda::ObsSpace &obsdb) {
  const std::vector<size_t> &obsCategories = obsdb.recnum();
  std::vector<size_t> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);
  splitter.groupBy(validObsCategories);
}

void TrackCheckUtils::groupObservationsByVariable(const Variable &variable,
                                                  const std::vector<size_t> &validObsIds,
                                                  RecursiveSplitter &splitter,
                                                  const ioda::ObsSpace &obsdb) {
  switch (obsdb.dtype(variable.group(), variable.variable())) {
  case ioda::ObsDtype::Integer:
    groupObservationsByTypedVariable<int>(variable, validObsIds, splitter, obsdb);
    break;

  case ioda::ObsDtype::String:
    groupObservationsByTypedVariable<std::string>(variable, validObsIds, splitter, obsdb);
    break;

  default:
    throw eckit::UserError(
          "Only integer and string variables may be used as station IDs", Here());
  }
}
template <typename VariableType>
void TrackCheckUtils::groupObservationsByTypedVariable(const Variable &variable,
                                                       const std::vector<size_t> &validObsIds,
                                                       RecursiveSplitter &splitter,
                                                       const ioda::ObsSpace &obsdb) {
  std::vector<VariableType> obsCategories(obsdb.nlocs());
  obsdb.get_db(variable.group(), variable.variable(), obsCategories);
  std::vector<VariableType> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);

  splitter.groupBy(validObsCategories);
}

void TrackCheckUtils::sortTracksChronologically(const std::vector<size_t> &validObsIds,
                                                RecursiveSplitter &splitter,
                                                const ioda::ObsSpace &obsdb) {
  std::vector<util::DateTime> times(obsdb.nlocs());
  obsdb.get_db("MetaData", "datetime", times);
  splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndexA, size_t obsIndexB)
  { return times[validObsIds[obsIndexA]] < times[validObsIds[obsIndexB]]; });
}

TrackCheckUtils::ObsGroupLocationTimes
  TrackCheckUtils::collectObservationsLocations(const ioda::ObsSpace &obsdb) {
  ObsGroupLocationTimes locationTimes;

  locationTimes.latitudes.resize(obsdb.nlocs());
  obsdb.get_db("MetaData", "latitude", locationTimes.latitudes);

  locationTimes.longitudes.resize(obsdb.nlocs());
  obsdb.get_db("MetaData", "longitude", locationTimes.longitudes);

  locationTimes.datetimes.resize(obsdb.nlocs());
  obsdb.get_db("MetaData", "datetime", locationTimes.datetimes);

  return locationTimes;
}

void TrackCheckUtils::flagRejectedObservations(const std::vector<bool> &isRejected,
                                               std::vector<std::vector<bool> >
                                               &flagged) {
  for (std::vector<bool> & variableFlagged : flagged)
    for (size_t obsId = 0; obsId < isRejected.size(); ++obsId)
      if (isRejected[obsId])
        variableFlagged[obsId] = true;
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
