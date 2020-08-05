/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/TrackCheckShip.h"

#include <Eigen/Core>  // for easy array addition/subtraction

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/filters/TrackCheckShipDiagnostics.h"
#include "ufo/filters/TrackCheckShipParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

/// \brief Estimate of speed as calculated in Ops_CheckShipTrack.
///
/// Speed is calculated between \p obs1 and \p obs2, accounting for spatial/temporal
/// uncertainty using the resolution values stored in \p options.
double TrackCheckShip::speedEstimate(
    const TrackCheckShip::TrackObservation &obs1,
    const TrackCheckShip::TrackObservation &obs2,
    const TrackCheckShipParameters& options) {
  util::Duration temporalDistance = abs(obs1.getTime() -
                                        obs2.getTime());
  util::Duration tempRes = options.temporalResolution;
  auto dist = distance(obs1, obs2);
  auto spatialRes = options.spatialResolution;
  double speedEst = 0.0;
  if (dist > spatialRes) {
    speedEst = (dist - 0.5 * spatialRes) /
        std::max(temporalDistance.toSeconds(),
                 (tempRes.toSeconds()));
  }
  return speedEst;
}

/// \brief Returns the angle in degrees (rounded to the nearest .1 degree)
/// formed by the track segment going from observation \p a to \p b to \p c.
float TrackCheckShip::angle(const TrackCheckShip::TrackObservation &a,
                                    const TrackCheckShip::TrackObservation &b,
                                    const TrackCheckShip::TrackObservation &c) {
  auto locA = a.getLocation();
  auto locB = b.getLocation();
  auto locC = c.getLocation();
  Eigen::Vector3f disp1{locB[0]-locA[0], locB[1]-locA[1], locB[2]-locA[2]};
  Eigen::Vector3f disp2{locC[0]-locB[0], locC[1]-locB[1], locC[2]-locB[2]};
  auto retValue = std::acos(disp1.dot(disp2) / (disp1.norm() * disp2.norm()));
  retValue *= static_cast<float>(Constants::rad2deg);
  return std::round(retValue * 10.0f) * 0.1f;
}

const TrackCheckShipDiagnostics* TrackCheckShip::diagnostics() const
{
  return diagnostics_.get();
}

TrackCheckShip::TrackObservation::TrackObservation(double latitude, double longitude,
                                                   const util::DateTime &time,
                                                   std::shared_ptr<TrackStatistics> const& ts)
  : obsLocationTime_(latitude, longitude, time), fullTrackStatistics_(ts) {}

TrackCheckShip::TrackCheckShip(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "TrackCheckShip: config = " << config << std::endl;
  options_.reset(new TrackCheckShipParameters());
  options_->deserialize(config);
  if (options_->testingMode)
    diagnostics_.reset(new TrackCheckShipDiagnostics());
}

// Required for the correct destruction of options_.
TrackCheckShip::~TrackCheckShip()
{}

void TrackCheckShip::print(std::ostream & os) const {
  os << "TrackCheckShip: config = " << config_ << std::endl;
}


/// \todo This section is not yet fully implemented. Current implementation includes separating
/// observations into tracks based on \p Station_Id, calculating distances, speeds, and
/// angles between observations, and incrementing track-specific counters should the
/// calculations produce unexpected values.
void TrackCheckShip::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  const std::vector<size_t> validObsIds = TrackCheckUtils::getValidObservationIds(apply, flags_);

  RecursiveSplitter splitter(validObsIds.size());
  TrackCheckUtils::groupObservationsByStation(validObsIds, splitter, config_, obsdb_);
  TrackCheckUtils::sortTracksChronologically(validObsIds, splitter, obsdb_);

  TrackCheckUtils::ObsGroupLocationTimes obsLocTime =
      TrackCheckUtils::collectObservationsLocations(obsdb_);
  for (auto track : splitter.multiElementGroups()) {
    std::vector<TrackObservation> trackObservations = collectTrackObservations(
          track.begin(), track.end(), validObsIds, obsLocTime);
    calculateTrackSegmentProperties(trackObservations, true);
  }
}

/// \returns a \p vector of \p TrackObservations that all hold a \p shared_ptr to an instance
/// of \p TrackStatistics, which holds all of the track-specific counters.
std::vector<TrackCheckShip::TrackObservation> TrackCheckShip::collectTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const TrackCheckUtils::ObsGroupLocationTimes &obsLocTime) const {
  std::vector<TrackObservation> trackObservations;
  trackObservations.reserve(trackObsIndicesEnd - trackObsIndicesBegin);
  std::shared_ptr<TrackStatistics> ts = std::make_shared<TrackStatistics>();
  for (std::vector<size_t>::const_iterator it = trackObsIndicesBegin;
       it != trackObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds[*it];
    trackObservations.push_back(TrackObservation(obsLocTime.latitudes[obsId],
                                                 obsLocTime.longitudes[obsId],
                                                 obsLocTime.datetimes[obsId], ts));
  }
  return trackObservations;
}

/// Calculates all of the statistics that require only two
/// adjacent \p TrackObservations,
/// storing within the righthand observation. This includes
/// distance between the two observations,
/// time difference between the observations, speed between the
/// observations, and if the
/// observations are recorded for the same time. Calls function to
/// increment track-wise counters based on
/// results.
void TrackCheckShip::TrackObservation::calculateTwoObservationValues(
    TrackObservation& prevObs, bool firstIteration,
    const TrackCheckShipParameters &options) {
  this->setDistance(TrackCheckShip::distance(prevObs, *this));
  (this->observationStatistics_.distance > options.spatialResolution) ?
        this->setSpeed(TrackCheckShip::speedEstimate(*this, prevObs, options)) :
        this->setSpeed(0.0);
  this->setTimeDifference(getTime() - prevObs.getTime());
  if (!(this->observationStatistics_.timeDifference.toSeconds())) {
    if (this->observationStatistics_.distance >
        options.spatialResolution) {
      prevObs.setSimultaneous(true);
    }
    setSimultaneous(true);
  }
  if (firstIteration) {
    adjustTwoObservationStatistics(options);
  }
}

/// Calculates all of the statistics that require three
/// adjacent \p TrackObservations,
/// storing within the middle observation. This includes
/// distance between two alternating observations,
/// speed between these alternating observations (if the middle
/// observation was not recorded), and the
/// angle formed by the three-observation track segment.
/// Calls function to increment track-wise counters based on results.
void TrackCheckShip::TrackObservation::calculateThreeObservationValues(
    const TrackObservation& prevObs, const TrackObservation& nextObs,
    bool firstIteration, const TrackCheckShipParameters &options) {
  this->setDistanceAveraged(TrackCheckShip::distance(prevObs, nextObs));
  this->setSpeedAveraged(speedEstimate(prevObs, nextObs, options));
  if (std::min(this->observationStatistics_.distance,
               nextObs.observationStatistics_.distance) >
      options.spatialResolution) {
    this->observationStatistics_.angle = angle(prevObs, *this, nextObs);
  }
  if (firstIteration) {
    adjustThreeObservationStatistics();
  }
}

/// This performs all of the necessary calculations based on the
/// observations' locations and timestamps by calling
/// \p calculateTwoObservationValues
/// and \p calculateThreeObservationValues for the
/// non-edge-case observations.
///
/// \todo After an observation rejection and removal,
/// certain calculations will need to be reperformed.
/// This feature is in progress.
void TrackCheckShip::calculateTrackSegmentProperties(
    std::vector<TrackObservation> &trackObservations,
    bool firstIteration) const {
  if (trackObservations.size()) {
    for (size_t obsIdx = 1; obsIdx < trackObservations.size(); obsIdx++) {
      TrackObservation &obs = trackObservations[obsIdx];
      TrackObservation &prevObs = trackObservations[obsIdx - 1];
      obs.calculateTwoObservationValues(prevObs, firstIteration, *options_);
      if (obsIdx > 1) {
        const TrackObservation &prevPrevObs = trackObservations[obsIdx - 2];
        prevObs.calculateThreeObservationValues(prevPrevObs, obs, firstIteration, *options_);
      }
      if (firstIteration && (obsIdx == trackObservations.size() - 1)) {
        int potentialDenominator = trackObservations.size() - 1 -
            obs.getFullTrackStatistics()->numShort_ - obs.getFullTrackStatistics()->numFast_;
        (obs.getFullTrackStatistics()->meanSpeed_) = (obs.getFullTrackStatistics()->sumSpeed_) /
            std::max(1, potentialDenominator);
      }
    }
    if (options_->testingMode) {
      std::vector<TrackCheckShip::ObservationStatistics> obsStats;
      for (size_t obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
        obsStats.push_back(trackObservations[obsIdx].getObservationStatistics());
      }
      auto trackStats = *(trackObservations[0].getFullTrackStatistics());
      diagnostics_->storeDiagnostics(std::make_pair(obsStats, trackStats));
    }
  }
}

/// \brief Keeps track of 0-distanced, short, and fast track segments,
/// as well as incrementing \p sumSpeed_ for normal track segments.
void TrackCheckShip::TrackObservation::adjustTwoObservationStatistics
(const TrackCheckShipParameters &options) const {
  util::Duration hour{"PT1H"};
  if (getObservationStatistics().distance == 0.0)
    getFullTrackStatistics()->numDist0_++;
  if (getObservationStatistics().timeDifference < hour) {
    getFullTrackStatistics()->numShort_++;
  } else if (getObservationStatistics().speed >= options.maxSpeed) {
    getFullTrackStatistics()->numFast_++;
  } else {
    getFullTrackStatistics()->sumSpeed_ += getObservationStatistics().speed;
  }
}

/// \brief Increments number of bends if angle is greater than or equal to 90 degrees
void TrackCheckShip::TrackObservation::adjustThreeObservationStatistics() const {
  if (getObservationStatistics().angle >= 90.0)
    getFullTrackStatistics()->numBends_++;
}

const TrackCheckShip::ObservationStatistics &
TrackCheckShip::TrackObservation::getObservationStatistics() const
{
  return observationStatistics_;
}

const std::shared_ptr<TrackCheckShip::TrackStatistics>
TrackCheckShip::TrackObservation::getFullTrackStatistics() const
{
  return fullTrackStatistics_;
}

/// Sets the calling observation as simultaneous and increments the track-wise
/// \p numSimultaneous_ counter.
void TrackCheckShip::TrackObservation::setSimultaneous(bool simul) {
  assert(simul);  // reverting a simultaneous obs to normal is not yet implemented
  if (!this->observationStatistics_.simultaneous)
    this->fullTrackStatistics_->numSimultaneous_++;
  this->observationStatistics_.simultaneous = simul;
}
void TrackCheckShip::TrackObservation::setDistance(double dist) {
  this->observationStatistics_.distance = dist;
}
void TrackCheckShip::TrackObservation::setTimeDifference(util::Duration tDiff) {
  this->observationStatistics_.timeDifference = tDiff;
}
void TrackCheckShip::TrackObservation::setSpeed(double speed) {
  this->observationStatistics_.speed = speed;
}
void TrackCheckShip::TrackObservation::setAngle(double angle) {
  this->observationStatistics_.angle = angle;
}
void TrackCheckShip::TrackObservation::setDistanceAveraged(double distAvg) {
  this->observationStatistics_.distanceAveraged = distAvg;
}
void TrackCheckShip::TrackObservation::setSpeedAveraged(double speedAvg) {
  this->observationStatistics_.speedAveraged = speedAvg;
}
}  // namespace ufo
