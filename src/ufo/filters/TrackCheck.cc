/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/TrackCheck.h"

#include <algorithm>
#include <cmath>
#include <map>
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
#include "ufo/filters/TrackCheckParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/utils/RecursiveSplitter.h"

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

enum Direction { FORWARD, BACKWARD, NUM_DIRECTIONS };

enum class CheckResult : char {
  FAILED = false,
  PASSED = true,
  SKIPPED
};

/// \brief Results of cross-checking an observation with another (a "buddy").
struct CheckResults {
  CheckResults() :
    isBuddyDistinct(false),
    speedCheckResult(CheckResult::SKIPPED),
    climbRateCheckResult(CheckResult::SKIPPED) {}

  bool isBuddyDistinct;
  CheckResult speedCheckResult;
  CheckResult climbRateCheckResult;
};

static const int NO_PREVIOUS_SWEEP = -1;

}  // namespace


/// \brief Locations of all observations processed by the track checking filter.
struct TrackCheck::ObsData {
  std::vector<float> latitudes;
  std::vector<float> longitudes;
  std::vector<util::DateTime> datetimes;
  std::vector<float> pressures;
};


/// \brief Attributes of an observation belonging to a track.
class TrackCheck::TrackObservation {
 public:
  TrackObservation(float latitude, float longitude, const util::DateTime &time, float pressure)
    : location_(pointFromLatLon(latitude, longitude)), time_(time), pressure_(pressure),
      rejectedInPreviousSweep_(false), rejectedBeforePreviousSweep_(false),
      numNeighborsVisitedInPreviousSweep_{NO_PREVIOUS_SWEEP, NO_PREVIOUS_SWEEP},
      numFailedChecks_(0), numChecks_(0)
  {}

  const Point &location() const { return location_; }
  const util::DateTime &time() const { return time_; }
  float pressure() const { return pressure_; }

  bool rejectedInPreviousSweep() const { return rejectedInPreviousSweep_; }
  bool rejectedBeforePreviousSweep() const { return rejectedBeforePreviousSweep_; }
  bool rejected() const { return rejectedInPreviousSweep_ || rejectedBeforePreviousSweep_; }

  float failedChecksFraction() const {
    return numChecks_ != 0 ? static_cast<float>(numFailedChecks_) / numChecks_ : 0.0f;
  }

  int numNeighborsVisitedInPreviousSweep(Direction dir) const {
    return numNeighborsVisitedInPreviousSweep_[dir];
  }

  void setNumNeighborsVisitedInPreviousSweep(Direction dir, int num) {
    numNeighborsVisitedInPreviousSweep_[dir] = num;
  }

  /// Estimates the instantaneous speed and climb rate by comparing this observation against
  /// \p buddyObs. Checks if these estimates are in the accepted ranges and if the two observations
  /// are far enough from each other to be considered "distinct".
  ///
  /// \param buddyObs Observation to compare against.
  /// \param options Track check options.
  /// \param maxValidSpeedAtPressure
  ///   Function mapping air pressure (in Pa) to the maximum realistic speed (in m/s).
  /// \param referencePressure
  ///   Pressure at which the maximum speed should be evaluated.
  ///
  /// \returns An object enapsulating the check results.
  CheckResults checkAgainstBuddy(const TrackObservation &buddyObs,
                                 const TrackCheckParameters &options,
                                 const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
                                 float referencePressure) const;

  void registerCheckResults(const CheckResults &result);
  void unregisterCheckResults(const CheckResults &result);

  void registerCheckResult(const CheckResult &result);
  void unregisterCheckResult(const CheckResult &result);

  void registerSweepOutcome(bool rejectedInSweep);

 private:
  Point location_;
  util::DateTime time_;
  float pressure_;
  bool rejectedInPreviousSweep_;
  bool rejectedBeforePreviousSweep_;
  int numNeighborsVisitedInPreviousSweep_[NUM_DIRECTIONS];
  int numFailedChecks_;
  int numChecks_;
};

CheckResults TrackCheck::TrackObservation::checkAgainstBuddy(
    const TrackObservation &buddyObs,
    const TrackCheckParameters &options,
    const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
    float referencePressure) const {
  CheckResults results;

  util::Duration temporalDistance = abs(buddyObs.time_ - time_);
  const float spatialDistance = distance(location_, buddyObs.location_);

  // Estimate the speed and check if it is within the allowed range
  const float conservativeSpeedEstimate =
      (spatialDistance - options.spatialResolution) /
      (temporalDistance + options.temporalResolution).toSeconds();
  const float maxSpeed = maxValidSpeedAtPressure(referencePressure);
  results.speedCheckResult = CheckResult(conservativeSpeedEstimate <= maxSpeed);

  // Estimate the climb rate and check if it is within the allowed range
  if (options.maxClimbRate.value() != boost::none) {
    const float pressureDiff = std::abs(pressure_ - buddyObs.pressure_);
    const float conservativeClimbRateEstimate =
        pressureDiff / (temporalDistance + options.temporalResolution).toSeconds();
    results.climbRateCheckResult =
        CheckResult(conservativeClimbRateEstimate <= *options.maxClimbRate.value());
  }

  const int resolutionMultiplier = options.distinctBuddyResolutionMultiplier;
  results.isBuddyDistinct =
      temporalDistance > resolutionMultiplier * options.temporalResolution &&
      spatialDistance > resolutionMultiplier * options.spatialResolution;

  return results;
}

void TrackCheck::TrackObservation::registerCheckResults(const CheckResults &results) {
  registerCheckResult(results.speedCheckResult);
  registerCheckResult(results.climbRateCheckResult);
  assert(numChecks_ >= 0);
  assert(numFailedChecks_ >= 0);
  assert(numFailedChecks_ <= numChecks_);
}

void TrackCheck::TrackObservation::unregisterCheckResults(const CheckResults &results) {
  unregisterCheckResult(results.speedCheckResult);
  unregisterCheckResult(results.climbRateCheckResult);
  assert(numChecks_ >= 0);
  assert(numFailedChecks_ >= 0);
  assert(numFailedChecks_ <= numChecks_);
}

void TrackCheck::TrackObservation::registerCheckResult(const CheckResult &result) {
  if (result != CheckResult::SKIPPED)
    ++numChecks_;
  if (result == CheckResult::FAILED)
    ++numFailedChecks_;
}

void TrackCheck::TrackObservation::unregisterCheckResult(const CheckResult &result) {
  if (result != CheckResult::SKIPPED)
    --numChecks_;
  if (result == CheckResult::FAILED)
    --numFailedChecks_;
}

void TrackCheck::TrackObservation::registerSweepOutcome(bool rejectedInSweep) {
  rejectedBeforePreviousSweep_ = rejected();
  rejectedInPreviousSweep_ = rejectedInSweep;
}


TrackCheck::TrackCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                       boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                       boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "TrackCheck: config = " << config_ << std::endl;

  options_.reset(new TrackCheckParameters());
  options_->deserialize(config);
}

// Required for the correct destruction of options_.
TrackCheck::~TrackCheck()
{}

void TrackCheck::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  const std::vector<size_t> validObsIds = getValidObservationIds(apply);

  RecursiveSplitter splitter(validObsIds.size());
  groupObservationsByStation(validObsIds, splitter);
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
    oops::Log::trace() << "TrackCheck: flagged? = " << flagged[0] << std::endl;
  }
}

std::vector<size_t> TrackCheck::getValidObservationIds(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && (*flags_)[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsId);
  return validObsIds;
}

void TrackCheck::groupObservationsByStation(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  if (options_->stationIdVariable.value() == boost::none) {
    if (obsdb_.obs_group_var().empty()) {
      // Observations were not grouped into records.
      // Assume all observations were taken during the same station.
      return;
    } else {
      groupObservationsByRecordNumber(validObsIds, splitter);
    }
  } else {
    groupObservationsByVariable(*options_->stationIdVariable.value(), validObsIds, splitter);
  }
}

void TrackCheck::groupObservationsByRecordNumber(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  const std::vector<size_t> &obsCategories = obsdb_.recnum();
  std::vector<size_t> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);
  splitter.groupBy(validObsCategories);
}

void TrackCheck::groupObservationsByVariable(
    const Variable &variable,
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  switch (obsdb_.dtype(variable.group(), variable.variable())) {
  case ioda::ObsDtype::Integer:
    groupObservationsByTypedVariable<int>(variable, validObsIds, splitter);
    break;

  case ioda::ObsDtype::String:
    groupObservationsByTypedVariable<std::string>(variable, validObsIds, splitter);
    break;

  default:
    throw eckit::UserError("Only integer and string variables may be used as station IDs", Here());
  }
}

template <typename VariableType>
void TrackCheck::groupObservationsByTypedVariable(
    const Variable &variable,
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  std::vector<VariableType> obsCategories(obsdb_.nlocs());
  obsdb_.get_db(variable.group(), variable.variable(), obsCategories);
  std::vector<VariableType> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);

  splitter.groupBy(validObsCategories);
}

void TrackCheck::sortTracksChronologically(const std::vector<size_t> &validObsIds,
                                           RecursiveSplitter &splitter) const {
  std::vector<util::DateTime> times(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "datetime", times);
  splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndexA, size_t obsIndexB)
  { return times[validObsIds[obsIndexA]] < times[validObsIds[obsIndexB]]; });
}

TrackCheck::ObsData TrackCheck::collectObsData() const {
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

PiecewiseLinearInterpolation TrackCheck::makeMaxSpeedByPressureInterpolation() const {
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

void TrackCheck::identifyRejectedObservationsInTrack(
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

std::vector<TrackCheck::TrackObservation> TrackCheck::collectTrackObservations(
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

TrackCheck::SweepResult TrackCheck::sweepOverObservations(
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
      const bool firstSweep = obs.numNeighborsVisitedInPreviousSweep(dir) == NO_PREVIOUS_SWEEP;
      const int numNeighborsVisitedInPreviousSweep =
          firstSweep ? 0 : obs.numNeighborsVisitedInPreviousSweep(dir);
      int numNewDistinctBuddiesToVisit = firstSweep ? options_->numDistinctBuddiesPerDirection : 0;

      auto getNthNeighbor = [&trackObservations, obsIdx, dir](int n) -> const TrackObservation* {
        const int neighborObsIdx = obsIdx + (dir == FORWARD ? n : -n);
        if (neighborObsIdx < 0 || neighborObsIdx >= trackObservations.size())
          return nullptr;  // We've reached the end of the track
        else
          return &trackObservations[neighborObsIdx];
      };

      float minPressureBetween = obs.pressure();
      int neighborIdx = 1;
      const TrackObservation *neighborObs = getNthNeighbor(neighborIdx);
      for (; neighborIdx <= numNeighborsVisitedInPreviousSweep && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        // Strictly speaking, neighborObs->pressure should be disregarded if neighborObs has already
        // been rejected. However, that would force us to check each pair of observations anew
        // whenever an observation between them is rejected, whereas as things stand, we only
        // need to "undo" checks against rejected observations.
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure());
        if (neighborObs->rejectedInPreviousSweep()) {
          CheckResults results = obs.checkAgainstBuddy(*neighborObs, *options_,
                                                       maxValidSpeedAtPressure, minPressureBetween);
          obs.unregisterCheckResults(results);
          if (results.isBuddyDistinct) {
            // The rejected distinct buddy needs to be replaced with another
            ++numNewDistinctBuddiesToVisit;
          }
        }
      }

      for (; numNewDistinctBuddiesToVisit > 0 && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure());
        if (!neighborObs->rejected()) {
          CheckResults results = obs.checkAgainstBuddy(*neighborObs, *options_,
                                                       maxValidSpeedAtPressure, minPressureBetween);
          obs.registerCheckResults(results);
          if (results.isBuddyDistinct)
            --numNewDistinctBuddiesToVisit;
        }
      }

      const int numNeighborsVisitedInThisSweep = neighborIdx - 1;
      assert(numNeighborsVisitedInThisSweep >= numNeighborsVisitedInPreviousSweep);
      obs.setNumNeighborsVisitedInPreviousSweep(dir, numNeighborsVisitedInThisSweep);
    }  // end of loop over directions

    failedChecksFraction[obsIdx] = obs.failedChecksFraction();
  }

  const float maxFailedChecksFraction = *std::max_element(failedChecksFraction.begin(),
                                                          failedChecksFraction.end());
  const float failedChecksThreshold = options_->rejectionThreshold * maxFailedChecksFraction;
  if (failedChecksThreshold <= 0)
    return SweepResult::NO_MORE_SWEEPS_REQUIRED;

  for (int obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
    const bool rejected = failedChecksFraction[obsIdx] > failedChecksThreshold;
    trackObservations[obsIdx].registerSweepOutcome(rejected);
  }

  return SweepResult::ANOTHER_SWEEP_REQUIRED;
}

void TrackCheck::flagRejectedTrackObservations(
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

void TrackCheck::flagRejectedObservations(const std::vector<bool> &isRejected,
                                          std::vector<std::vector<bool> > &flagged) const {
  for (std::vector<bool> & variableFlagged : flagged)
    for (size_t obsId = 0; obsId < isRejected.size(); ++obsId)
      if (isRejected[obsId])
        variableFlagged[obsId] = true;
}

void TrackCheck::print(std::ostream & os) const {
  os << "TrackCheck: config = " << config_ << std::endl;
}

}  // namespace ufo
