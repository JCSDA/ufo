/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/AircraftTrackCheck.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/Sphere.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/sqr.h"
#include "ufo/filters/AircraftTrackCheckParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace util {
inline Duration abs(const Duration &duration) {
  return duration.toSeconds() >= 0 ? duration : -duration;
}
}  // namespace util

namespace ufo {

namespace {

const int worldDim = 3;
typedef std::array<float, worldDim> Point;

float distance2(const Point &a, const Point &b) {
  float sum = 0;
  for (size_t i = 0; i < a.size(); ++i)
    sum += util::sqr(a[i] - b[i]);
  return sum;
}

float distance(const Point &a, const Point &b) {
  return std::sqrt(distance2(a, b));
}

Point pointFromLatLon(float latitude, float longitude) {
  // This local copy is needed because convertSphericalToCartesian takes the first parameter by
  // reference, but Constants::mean_earth_rad has no out-of-line definition.
  const double meanEarthRadius = Constants::mean_earth_rad;
  eckit::geometry::Point3 eckitPoint;
  eckit::geometry::Sphere::convertSphericalToCartesian(
        meanEarthRadius, eckit::geometry::Point2(longitude, latitude), eckitPoint);
  Point Point;
  std::copy(eckitPoint.begin(), eckitPoint.end(), Point.begin());
  return Point;
}

const int UNKNOWN = -1;
enum Direction { FORWARD, BACKWARD, NUM_DIRECTIONS };

}  // namespace

/// \brief Locations of all observations processed by the track checking filter.
struct AircraftTrackCheck::ObsData {
  std::vector<float> latitudes;
  std::vector<float> longitudes;
  std::vector<util::DateTime> datetimes;
  std::vector<float> pressures;
};

/// \brief Attributes of an observation belonging to a track.
struct AircraftTrackCheck::TrackObservation {
  TrackObservation(float latitude, float longitude, const util::DateTime &time_, float pressure_)
    : location(pointFromLatLon(latitude, longitude)), time(time_), pressure(pressure_),
      rejectedInPreviousSweep(false), rejectedBeforePreviousSweep(false),
      numNeighborsToVisit{UNKNOWN, UNKNOWN},
      numFailedChecks(0), numChecks(0)
  {}

  bool rejected() const { return rejectedInPreviousSweep || rejectedBeforePreviousSweep; }

  Point location;
  util::DateTime time;
  float pressure;
  bool rejectedInPreviousSweep;
  bool rejectedBeforePreviousSweep;
  int numNeighborsToVisit[NUM_DIRECTIONS];
  int numFailedChecks;
  int numChecks;
};

/// \brief Results of cross-checking an observation with another (a "buddy").
struct AircraftTrackCheck::CheckResult {
  bool isBuddyDistinct;
  bool speedCheckPassed;
  bool climbRateCheckPassed;
};

AircraftTrackCheck::AircraftTrackCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                       boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                       boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "AircraftTrackCheck: config = " << config_ << std::endl;

  options_.reset(new AircraftTrackCheckParameters());
  options_->deserialize(config);
}

// Required for the correct destruction of options_.
AircraftTrackCheck::~AircraftTrackCheck()
{}

void AircraftTrackCheck::applyFilter(const std::vector<bool> & apply,
                                     const Variables & filtervars,
                                     std::vector<std::vector<bool>> & flagged) const {
  const std::vector<size_t> validObsIds = getValidObservationIds(apply);

  RecursiveSplitter splitter(validObsIds.size());
  groupObservationsByFlightId(validObsIds, splitter);
  sortTracksChronologically(validObsIds, splitter);

  ObsData obsData = collectObsData();
  PiecewiseLinearInterpolation maxSpeedByPressure = makeMaxSpeedByPressureInterpolation();

  std::vector<bool> isRejected(apply.size(), false);
  for (auto track : splitter.multiElementGroups()) {
    identifyRejectedObservationsInTrack(track.begin(), track.end(), validObsIds,
                                        obsData, maxSpeedByPressure, isRejected);
  }
  flagRejectedObservations(isRejected, flagged);

  if (filtervars.size() != 0) {
    oops::Log::trace() << "AircraftTrackCheck: flagged? = " << flagged[0] << std::endl;
  }
}

std::vector<size_t> AircraftTrackCheck::getValidObservationIds(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && (*flags_)[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsId);
  return validObsIds;
}

void AircraftTrackCheck::groupObservationsByFlightId(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  ioda::ObsDataVector<int> obsDataVector(obsdb_, options_->flightIdVariable.value().variable(),
                                         options_->flightIdVariable.value().group());
  const auto &flightId = obsDataVector[0];

  std::vector<int> validObsCategories(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex)
    validObsCategories[validObsIndex] = flightId[validObsIds[validObsIndex]];
  splitter.groupBy(validObsCategories);
}

void AircraftTrackCheck::sortTracksChronologically(const std::vector<size_t> &validObsIds,
                                                   RecursiveSplitter &splitter) const {
  std::vector<util::DateTime> times(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "datetime", times);
  splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndexA, size_t obsIndexB)
                        { return times[validObsIds[obsIndexA]] < times[validObsIds[obsIndexB]]; });
}

AircraftTrackCheck::ObsData AircraftTrackCheck::collectObsData() const {
  ObsData obsData;

  obsData.latitudes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "latitude", obsData.latitudes);

  obsData.longitudes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "longitude", obsData.longitudes);

  obsData.datetimes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "datetime", obsData.datetimes);

  obsData.pressures.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "air_pressure", obsData.pressures);

  return obsData;
}

PiecewiseLinearInterpolation AircraftTrackCheck::makeMaxSpeedByPressureInterpolation() const {
  const std::map<float, float> &maxSpeedInterpolationPoints =
      options_->maxSpeedInterpolationPoints.value();

  std::vector<double> pressures, maxSpeeds;
  pressures.reserve(maxSpeedInterpolationPoints.size());
  maxSpeeds.reserve(maxSpeedInterpolationPoints.size());

  for (const auto &pressureAndMaxSpeed : maxSpeedInterpolationPoints) {
    pressures.push_back(pressureAndMaxSpeed.first);
    // The interpolator needs to produce speeds in km/s rather than m/s because observation
    // locations are expressed in kilometers.
    const int metersPerKm = 1000;
    maxSpeeds.push_back(pressureAndMaxSpeed.second / metersPerKm);
  }

  return PiecewiseLinearInterpolation(std::move(pressures), std::move(maxSpeeds));
}

void AircraftTrackCheck::identifyRejectedObservationsInTrack(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const ObsData &obsData,
    const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
    std::vector<bool> &isRejected) const {

  std::vector<TrackObservation> trackObservations = collectTrackObservations(
        trackObsIndicesBegin, trackObsIndicesEnd, validObsIds, obsData);
  std::vector<float> workspace;

  while (sweepOverObservations(trackObservations, maxValidSpeedAtPressure, workspace) ==
         SweepResult::ANOTHER_SWEEP_REQUIRED) {
    // can't exit the loop yet
  }

  flagRejectedTrackObservations(trackObsIndicesBegin, trackObsIndicesEnd,
                                validObsIds, trackObservations, isRejected);
}

std::vector<AircraftTrackCheck::TrackObservation> AircraftTrackCheck::collectTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const ObsData &obsData) const {
  std::vector<TrackObservation> trackObservations;
  trackObservations.reserve(trackObsIndicesEnd - trackObsIndicesBegin);
  for (std::vector<size_t>::const_iterator it = trackObsIndicesBegin;
       it != trackObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds[*it];
    trackObservations.push_back(TrackObservation(obsData.latitudes[obsId],
                                                 obsData.longitudes[obsId],
                                                 obsData.datetimes[obsId],
                                                 obsData.pressures[obsId]));
  }
  return trackObservations;
}

AircraftTrackCheck::SweepResult AircraftTrackCheck::sweepOverObservations(
    std::vector<TrackObservation> &trackObservations,
    const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
    std::vector<float> &workspace) const {

  std::vector<float> &failedChecksFraction = workspace;
  failedChecksFraction.assign(trackObservations.size(), 0.0f);

  for (int obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
    TrackObservation &obs = trackObservations[obsIdx];
    if (obs.rejected())
      continue;

    for (Direction dir : { FORWARD, BACKWARD}) {
      const bool firstSweep = obs.numNeighborsToVisit[dir] == UNKNOWN;      
      const int numNeighborsVisitedInPreviousSweep = firstSweep ? 0 : obs.numNeighborsToVisit[dir];
      int numNewDistinctBuddiesToVisit = firstSweep ? options_->numDistinctBuddiesPerDirection : 0;

      auto getNthNeighbor = [&trackObservations, obsIdx, dir](int n) -> const TrackObservation* {
        const int neighborObsIdx = obsIdx + (dir == FORWARD ? n : -n);
        if (neighborObsIdx < 0 || neighborObsIdx >= trackObservations.size())
          return nullptr; // We've reached the end of the track
        else
          return &trackObservations[neighborObsIdx];
      };

      float minPressureBetween = obs.pressure;
      int neighborIdx = 1;
      const TrackObservation *neighborObs = getNthNeighbor(neighborIdx);
      for (; neighborIdx <= numNeighborsVisitedInPreviousSweep && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        // Strictly speaking, neighborObs->pressure should be disregarded if neighborObs has already
        // been rejected. However, that would force us to check each pair of observations anew
        // whenever an observation lying in-between is rejected, whereas as things stand, we only
        // need to revisit
        // sandwiching a rejected observation from
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure);
        if (neighborObs->rejectedInPreviousSweep) {
          CheckResult result = checkObservationPair(obs, *neighborObs,
                                                    maxValidSpeedAtPressure, minPressureBetween);
          unregisterCheckResult(obs, result);
          if (result.isBuddyDistinct) {
            // We must replace the rejected distinct buddy with another
            ++numNewDistinctBuddiesToVisit;
          }
        }
      }

      for (; numNewDistinctBuddiesToVisit > 0 && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure);
        if (!neighborObs->rejected()) {
          CheckResult result = checkObservationPair(obs, *neighborObs,
                                                    maxValidSpeedAtPressure, minPressureBetween);
          registerCheckResult(obs, result);
          if (result.isBuddyDistinct)
            --numNewDistinctBuddiesToVisit;
        }
      }

      obs.numNeighborsToVisit[dir] = neighborIdx - 1;
      assert(obs.numNeighborsToVisit[dir] >= numNeighborsVisitedInPreviousSweep);
    }  // end of loop over directions

    failedChecksFraction[obsIdx] =
        obs.numChecks != 0 ? static_cast<float>(obs.numFailedChecks) / obs.numChecks : 0;
  }

  const float maxFailedChecksFraction = *std::max_element(failedChecksFraction.begin(),
                                                          failedChecksFraction.end());
  const float failedChecksThreshold = options_->rejectionThreshold * maxFailedChecksFraction;
  if (failedChecksThreshold <= 0)
    return SweepResult::NO_MORE_SWEEPS_REQUIRED;

  for (int obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
    TrackObservation &obs = trackObservations[obsIdx];
    obs.rejectedBeforePreviousSweep = obs.rejected();
    obs.rejectedInPreviousSweep = failedChecksFraction[obsIdx] > failedChecksThreshold;
  }

  return SweepResult::ANOTHER_SWEEP_REQUIRED;
}

AircraftTrackCheck::CheckResult AircraftTrackCheck::checkObservationPair(
    const TrackObservation &obs, const TrackObservation &buddyObs,
    const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
    float minPressureBetween) const {
  CheckResult result;

  util::Duration temporalDistance = abs(buddyObs.time - obs.time);
  const float spatialDistance = distance(obs.location, buddyObs.location);

  // Estimate the speed and check if it is within the allowed range
  const float conservativeSpeedEstimate =
      (spatialDistance - options_->spatialResolution) /
      (temporalDistance + options_->temporalResolution).toSeconds();
  const float maxSpeed = maxValidSpeedAtPressure(minPressureBetween);
  result.speedCheckPassed = (conservativeSpeedEstimate <= maxSpeed);

  // Estimate the climb rate and check if it is within the allowed range
  const float pressureDiff = std::abs(obs.pressure - buddyObs.pressure);
  const float conservativeClimbRateEstimate =
      pressureDiff / (temporalDistance + options_->temporalResolution).toSeconds();
  result.climbRateCheckPassed = (conservativeClimbRateEstimate <= options_->maxClimbRate);

  const int resolutionMultiplier = options_->distinctBuddyResolutionMultiplier;
  result.isBuddyDistinct =
      temporalDistance > resolutionMultiplier * options_->temporalResolution &&
      spatialDistance > resolutionMultiplier * options_->spatialResolution;

  return result;
}

void AircraftTrackCheck::registerCheckResult(TrackObservation &obs,
                                               const CheckResult &result) const {
  obs.numChecks += 2;
  if (!result.speedCheckPassed)
    ++obs.numFailedChecks;
  if (!result.climbRateCheckPassed)
    ++obs.numFailedChecks;
  assert(obs.numChecks >= 0);
  assert(obs.numFailedChecks >= 0);
  assert(obs.numFailedChecks <= obs.numChecks);
}

void AircraftTrackCheck::unregisterCheckResult(TrackObservation &obs,
                                               const CheckResult &result) const {
  obs.numChecks -= 2;
  if (!result.speedCheckPassed)
    --obs.numFailedChecks;
  if (!result.climbRateCheckPassed)
    --obs.numFailedChecks;
  assert(obs.numChecks >= 0);
  assert(obs.numFailedChecks >= 0);
  assert(obs.numFailedChecks <= obs.numChecks);
}

void AircraftTrackCheck::flagRejectedTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const std::vector<TrackObservation> &trackObservations,
    std::vector<bool> &isRejected) const {
  auto trackObsIndexIt = trackObsIndicesBegin;
  auto trackObsIt = trackObservations.begin();
  for (; trackObsIndexIt != trackObsIndicesEnd; ++trackObsIndexIt, ++trackObsIt)
    if (trackObsIt->rejected())
      isRejected[validObsIds[*trackObsIndexIt]] = true;
}

void AircraftTrackCheck::flagRejectedObservations(const std::vector<bool> &isRejected,
                                                  std::vector<std::vector<bool> > &flagged) const {
  for (std::vector<bool> & variableFlagged : flagged)
    for (size_t obsId = 0; obsId < isRejected.size(); ++obsId)
      if (isRejected[obsId])
        variableFlagged[obsId] = true;
}

void AircraftTrackCheck::print(std::ostream & os) const {
  os << "AircraftTrackCheck: config = " << config_ << std::endl;
}

}  // namespace ufo
