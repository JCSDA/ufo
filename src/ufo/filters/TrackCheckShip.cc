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
#include <functional>
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
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/TrackCheckShipDiagnostics.h"
#include "ufo/filters/TrackCheckShipParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

TrackCheckShip::TrackCheckShip(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "TrackCheckShip: config = " << config << std::endl;
  options_.reset(new TrackCheckShipParameters());
  options_->deserialize(config);
  assert(options_->maxSpeed.value() > 0);
  if (options_->testingMode)
    diagnostics_.reset(new TrackCheckShipDiagnostics());
}

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
    speedEst *= 1000.0;  // convert units from km/s to m/s
  }
  return speedEst;
}

/// \brief Returns the angle in degrees (rounded to the nearest .1 degree)
/// between displacement vectors going from observations \p a to \p b
/// and \p b to \p c.
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

TrackCheckShip::TrackObservation::TrackObservation(
    double latitude, double longitude,
    const util::DateTime &time,
    std::shared_ptr<TrackStatistics> const& ts,
    std::shared_ptr<TrackCheckUtils::CheckCounter> const &checkCounter,
    size_t observationNumber)
  : obsLocationTime_(latitude, longitude, time), fullTrackStatistics_(ts),
    observationNumber_(observationNumber), rejected_(false) {}

size_t TrackCheckShip::TrackObservation::getObservationNumber() const
{
  return observationNumber_;
}

// Required for the correct destruction of options_.
TrackCheckShip::~TrackCheckShip()
{}

/// Detects whether each observation in a track has been rejected, and marks
/// the corresponding space in \p isRejected as true if so.
/// \param trackObsIndicesBegin the begin iterator of the track
/// \param trackObsIndicesEnd the end iterator
/// \param validObsIds the vector marking all of the observations' global
/// positions
/// \param trackObservations the full vector of observations within
/// the single track
/// \param isRejected the boolean vector whose indices correspond to all of
/// the rejected observations in the full input dataset
void TrackCheckShip::flagRejectedTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const std::vector<TrackObservation> &trackObservations,
    std::vector<bool> &isRejected) const {
  auto trackObsIndexIt = trackObsIndicesBegin;
  auto trackObsIt = trackObservations.begin();
  for (; trackObsIndexIt != trackObsIndicesEnd; ++trackObsIndexIt, ++trackObsIt)
    isRejected[validObsIds[*trackObsIndexIt]] = trackObsIt->rejected();
}

void TrackCheckShip::print(std::ostream & os) const {
  os << "TrackCheckShip: config = " << config_ << std::endl;
}


/// The filter begins with separating observations into tracks based on \p Station_Id. Various
/// calculations are then performed between consecutive observations: distances/speeds between
/// observation pairs, and angles between the displacement vectors formed by triplets of
/// consecutive observations. Track-specific counters are incremented for various metrics that
/// indicate an unexpected track path.
///
/// Based on these counter values, the filter may ignore the current track and move on to the next
/// one. Using the described calculations and additional ones if necessary, the observation
/// directly before and/or after the fastest segment is chosen for rejection. This will continue
/// until all remaining segments that link accepted observations are either:
/// 1. slower than a user-defined "max speed (m/s)" where the angles between the fastest accepted
/// segment's displacement vector and its two neighbouring segments' displacement vectors are
/// both less than 90 degrees.
/// 2. slower than 80 percent of this user-defined speed.
///
/// The full track will be rejected if the number of observations removed is more than the
/// user-defined "rejection threshold".

void TrackCheckShip::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(options_->stationIdVariable, obsdb_);

  const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply, *flags_);

  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  TrackCheckUtils::sortTracksChronologically(validObsIds, obsAccessor, splitter);

  TrackCheckUtils::ObsGroupLocationTimes obsLocTime =
      TrackCheckUtils::collectObservationsLocations(obsAccessor);

  std::vector<bool> isRejected(obsLocTime.latitudes.size(), false);
  size_t trackNumber = 0;
  for (auto track : splitter.multiElementGroups()) {
    trackNumber++;
    std::string stationId = std::to_string(trackNumber);
    std::vector<TrackObservation> trackObservations = collectTrackObservations(
          track.begin(), track.end(), validObsIds, obsLocTime);
    std::vector<std::reference_wrapper<TrackObservation>> trackObservationsReferences;
    trackObservationsReferences.reserve(trackObservations.size());
    std::transform(trackObservations.begin(), trackObservations.end(),
                   std::back_inserter(trackObservationsReferences),
                   [](TrackObservation &obs) {
      return std::ref<TrackObservation>(obs); });
    calculateTrackSegmentProperties(trackObservationsReferences,
                                    CalculationMethod::FIRSTITERATION);
    if (!trackObservationsReferences.empty() &&
        this->options_->earlyBreakCheck &&
        TrackCheckShip::earlyBreak(trackObservationsReferences, stationId)) {
      continue;
    }
      bool firstIterativeRemoval = true;
      while (trackObservationsReferences.size() >= 3) {
        // Initial loop: fastest (as determined by set of comparisons) observation removed
        // until all segments show slower speed than max threshold
        auto maxSpeedReferenceIterator = std::max_element(
              trackObservationsReferences.begin(), trackObservationsReferences.end(),
              [](TrackObservation a, TrackObservation b) {
            return a.getObservationStatistics().speed <
            b.getObservationStatistics().speed;});
        auto maxSpeedValue = maxSpeedReferenceIterator->get().getObservationStatistics().speed;
        if (maxSpeedValue <= (0.8 * options_->maxSpeed.value())) {
          break;
        } else if (maxSpeedValue < options_->maxSpeed.value()) {
          auto maxSpeedAngle = std::max(
                (maxSpeedReferenceIterator - 1)->get().getObservationStatistics().angle,
                maxSpeedReferenceIterator->get().getObservationStatistics().angle);
          if (maxSpeedAngle <= 90.0) {
            break;
          }
        }
        removeFaultyObservation(
              trackObservationsReferences, maxSpeedReferenceIterator, firstIterativeRemoval,
              stationId);
        firstIterativeRemoval = false;
        calculateTrackSegmentProperties(trackObservationsReferences, CalculationMethod::MAINLOOP);
      }
      auto rejectedCount = std::count_if(trackObservations.begin(), trackObservations.end(),
                    [](const TrackObservation& a) {return a.rejected();});
      if (rejectedCount >= options_->rejectionThreshold.value() * trackObservations.size()) {
        oops::Log::trace() << "CheckShipTrack: track " << stationId << " NumRej " <<
                              rejectedCount << " out of " << trackObservations.size() <<
                              " reports rejected. *** Reject whole track ***\n";
        for (TrackObservation &obs : trackObservations)
          obs.setRejected(true);
      }
      flagRejectedTrackObservations(track.begin(), track.end(),
                                    validObsIds, trackObservations, isRejected);
  }
  obsAccessor.flagRejectedObservations(isRejected, flagged);
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
  std::shared_ptr<TrackStatistics> trackStatistics(new TrackStatistics());
  std::shared_ptr<TrackCheckUtils::CheckCounter> checkCounter(new TrackCheckUtils::CheckCounter);
  size_t observationNumber = 0;
  for (std::vector<size_t>::const_iterator it = trackObsIndicesBegin;
       it != trackObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds[*it];
    trackObservations.push_back(TrackObservation(obsLocTime.latitudes[obsId],
                                                 obsLocTime.longitudes[obsId],
                                                 obsLocTime.datetimes[obsId], trackStatistics,
                                                 checkCounter,
                                                 observationNumber));
    observationNumber++;
  }
  return trackObservations;
}

/// \brief \returns true if at least half of the track segments have
/// incremented the relevant rejection counters
///
/// This filter is best for mostly-accurate observation tracks. If the track has
/// many "errors" (as indicated by the counters that are incremented before
/// any observations are removed), it will stop early by default.
///
/// Sometimes caused by two ships with same callsign or various reports with wrong time,
/// the check gives up. This is particularly a problem with WOD01 data - case studies
/// suggest that most suspect data is reasonable.
bool TrackCheckShip::earlyBreak(const std::vector<std::reference_wrapper<TrackObservation>>
                                &trackObs, const std::string trackId) const {
  bool breakResult = false;
  const auto& trackStats = *(trackObs[0].get().getFullTrackStatistics());
  // if at least half of the track segments have a time difference of less than an hour
  // (if non-buoy), are faster than a configured maximum speed, or exhibit at least a 90
  // degree bend
  if ((2 * ((options_->inputCategory.value() != 1)  // 1 is input category of buoy
            * trackStats.numShort_ + trackStats.numFast_) + trackStats.numBends_)
      >= (trackObs.size() - 1)) {
    oops::Log::trace() << "ShipTrackCheck: " << trackId << "\n" <<
                          "Time difference < 1 hour: " << trackStats.numShort_ << "\n" <<
                          "Fast: " << trackStats.numFast_ << "\n" <<
                          "Bends: " << trackStats.numBends_ << "\n" <<
                          "Total observations: " << trackObs.size() << "\n" <<
                          "Track was not checked." << std::endl;

    breakResult = true;
  }
  if (options_->testingMode)
    diagnostics_->storeEarlyBreakResult(breakResult);
  return breakResult;
}

/// \brief Chooses which of the observations surrounding the speediest segment to remove,
/// flagging it accordingly.
void TrackCheckShip::removeFaultyObservation(
    std::vector<std::reference_wrapper<TrackObservation>> &track,
    const std::vector<std::reference_wrapper<TrackObservation>>::iterator
    &observationAfterFastestSegment,
    bool firstIterativeRemoval, const std::string trackId) const {
  int errorCategory = 0;
  util::Duration four_days{"P4D"};
  auto rejectedObservation = observationAfterFastestSegment;
  // lambda function to "fail" an observation that should be rejected
  auto fail = [&rejectedObservation](
      const std::vector<std::reference_wrapper<TrackObservation>>::iterator &iter) {
    iter->get().setRejected(true);
    rejectedObservation = iter;
  };
  auto neighborObservationStatistics = [&observationAfterFastestSegment](int index) ->
      const ObservationStatistics & {
    return (observationAfterFastestSegment + index)->get().getObservationStatistics();
  };
  auto meanSpeed = observationAfterFastestSegment->get().getFullTrackStatistics()->meanSpeed_;
  if (observationAfterFastestSegment == track.begin() + 1) {
    // Decide whether ob 0 or 1 agrees best with ob 2
    if (neighborObservationStatistics(0).speedAveraged <=
        options_->maxSpeed &&
        (neighborObservationStatistics(1).speed >
         options_->maxSpeed ||
         neighborObservationStatistics(1).angle > 45.0)) {
      fail(observationAfterFastestSegment);
      errorCategory = 2;
    } else {
      fail(observationAfterFastestSegment - 1);
      errorCategory = 1;
    }
  } else if (observationAfterFastestSegment == track.end() - 1) {
    if (neighborObservationStatistics(-1).speedAveraged <=
        options_->maxSpeed &&
        (neighborObservationStatistics(-1).speed >
         options_->maxSpeed ||
         neighborObservationStatistics(-2).angle > 45.0)) {
      fail(observationAfterFastestSegment - 1);
      errorCategory = 2;
    } else {
      fail(observationAfterFastestSegment);
      errorCategory = 1;
    }
  } else if (neighborObservationStatistics(-1).speed >
             options_->maxSpeed) {
    fail(observationAfterFastestSegment - 1);
    errorCategory = 4;
    //  Category 4: both segments surrounding observation have excessive speed
  } else if (neighborObservationStatistics(1).speed >
             options_->maxSpeed) {
    fail(observationAfterFastestSegment);
    errorCategory = 4;
  } else if (neighborObservationStatistics(0).speedAveraged >
             options_->maxSpeed) {
    fail(observationAfterFastestSegment - 1);
    errorCategory = 5;
    // Category 5: observation before fastest segment would still begin a fast segment
    // if observation after fastest segment were removed
  } else if (neighborObservationStatistics(-1).speedAveraged >
             options_->maxSpeed) {
    fail(observationAfterFastestSegment);
    errorCategory = 5;
    // Category 5: observation after fastest segment would still end a fast segment if
    // observation before fastest segment were removed
  } else if (neighborObservationStatistics(-1).angle >
             (45.0 + neighborObservationStatistics(0).angle)) {
    fail(observationAfterFastestSegment - 1);
    errorCategory = 6;
  } else if (neighborObservationStatistics(0).angle >
             45.0 + neighborObservationStatistics(-1).angle) {
    fail(observationAfterFastestSegment);
    errorCategory = 6;
  } else if (neighborObservationStatistics(-2).angle > 45.0 &&
             neighborObservationStatistics(-2).angle >
             neighborObservationStatistics(1).angle) {
    fail(observationAfterFastestSegment - 1);
    errorCategory = 7;
  } else if (neighborObservationStatistics(1).angle > 45.0) {
    fail(observationAfterFastestSegment);
    errorCategory = 7;
  } else if (neighborObservationStatistics(-1).speed <
             0.5 * std::min(
               neighborObservationStatistics(1).speed,
               meanSpeed)) {
    fail(observationAfterFastestSegment - 1);
    errorCategory = 8;
  } else if (neighborObservationStatistics(1).speed <
             0.5 * std::min(neighborObservationStatistics(-1).speed, meanSpeed)) {
    fail(observationAfterFastestSegment);
    errorCategory = 8;
  } else {
    double distanceSum = 0.0;
    distanceSum = std::accumulate(observationAfterFastestSegment - 1,
                    observationAfterFastestSegment + 2, distanceSum,
                    [](double summed, TrackObservation& obs) {
      return summed + obs.getObservationStatistics().distance;
    });

    double distancePrevObsOmitted =
          neighborObservationStatistics(-1).distanceAveraged +
        neighborObservationStatistics(1).distance;
    double distanceCurrentObsOmitted =
          neighborObservationStatistics(-1).distance +
        neighborObservationStatistics(0).distanceAveraged;
    util::Duration timeSum = (observationAfterFastestSegment + 1)->get().getTime() - (
          observationAfterFastestSegment - 2)->get().getTime();
    if (options_->testingMode.value()) {
      diagnostics_->storeDistanceSum(distanceSum);
      diagnostics_->storeDistancePrevObsOmitted(distancePrevObsOmitted);
      diagnostics_->storeDistanceCurrentObsOmitted(distanceCurrentObsOmitted);
      double timeDouble = timeSum.toSeconds();
      diagnostics_->storeTimeSum(timeDouble);
    }
    if (distancePrevObsOmitted < distanceCurrentObsOmitted - std::max(
          options_->spatialResolution.value(), 0.1 * distanceSum)) {
      fail(observationAfterFastestSegment - 1);
      errorCategory = 9;
    } else if (distanceCurrentObsOmitted < (
                 distancePrevObsOmitted - std::max(
                   options_->spatialResolution.value(), 0.1 * distanceSum))) {
      fail(observationAfterFastestSegment);
      errorCategory = 9;
    } else if (timeSum <= four_days && timeSum.toSeconds() > 0 &&
               std::min(distancePrevObsOmitted, distanceCurrentObsOmitted) > 0.0) {
      double previousSegmentDistanceProportion =
          // Prev segment dist/(prev segment distance + distAveragedCurrentObservation)
          neighborObservationStatistics(-1).distance /
          distanceCurrentObsOmitted;
      double previousObservationDistanceAveragedProportion =
          // Previous observation distAveraged / (prev obs distAveraged + next segment distance)
          neighborObservationStatistics(-1).
          distanceAveraged / distancePrevObsOmitted;
      double previousSegmentTimeProportion =
          static_cast<double>(((observationAfterFastestSegment - 1)->get().getTime() -
           (observationAfterFastestSegment - 2)->get().getTime()).toSeconds()) /
          timeSum.toSeconds();
      double previousAndFastestSegmentTimeProportion =
          static_cast<double>(((observationAfterFastestSegment->get().getTime()) -
            (observationAfterFastestSegment - 2)->get().getTime()).toSeconds()) /
          timeSum.toSeconds();
      if (options_->testingMode.value()) {
        diagnostics_->storePreviousSegmentDistanceProportion(previousSegmentDistanceProportion);
        diagnostics_->storePreviousObservationDistanceAveragedProportion(
              previousObservationDistanceAveragedProportion);
        diagnostics_->storePreviousSegmentTimeProportion(previousSegmentTimeProportion);
        diagnostics_->storePreviousAndFastestSegmentTimeProportion(
              previousAndFastestSegmentTimeProportion);
      }
      if (std::abs(previousSegmentDistanceProportion - previousSegmentTimeProportion) >
          0.1 + std::abs(previousObservationDistanceAveragedProportion -
                         previousAndFastestSegmentTimeProportion)) {
        // previous segment's spatial and temporal lengths are significantly more disproportionate
        // when compared to a larger portion of track
        // than the equivalents for the post-fastest segment
        fail(observationAfterFastestSegment - 1);
        errorCategory = 10;
      } else if (std::abs(previousObservationDistanceAveragedProportion -
                          previousAndFastestSegmentTimeProportion) > 0.1 +
                 std::abs(previousSegmentDistanceProportion - previousSegmentTimeProportion)) {
        // next segment after fastest has spatial and temporal lengths that are significantly more
        // disproportionate than the segment before the fastest
        fail(observationAfterFastestSegment);
        errorCategory = 10;
      } else {
        fail(observationAfterFastestSegment);
      }
      oops::Log::trace() << "CheckShipTrack: proportions " << previousSegmentDistanceProportion <<
                            " " << previousSegmentTimeProportion <<
                            " " << previousObservationDistanceAveragedProportion << " "
                         << previousAndFastestSegmentTimeProportion << " speeds: " << meanSpeed
                         << " " << neighborObservationStatistics(-1).speed << " " <<
                            neighborObservationStatistics(0).speed << " " <<
                            neighborObservationStatistics(1).speed << " [m/s]" << std::endl;
    }
    if (errorCategory == 9 || std::min(distancePrevObsOmitted, distanceCurrentObsOmitted) == 0.0) {
      oops::Log::trace() << "CheckShipTrack: Dist check, station id: " <<
                            trackId << std::endl <<
                            " error category: " << errorCategory << std::endl <<
                            " distances: " << distanceSum * 0.001 << " " <<
                            distancePrevObsOmitted * 0.001 << " " <<
                            distanceCurrentObsOmitted * 0.001 << " " <<
                            (distancePrevObsOmitted - distanceCurrentObsOmitted) * 0.001 <<
                            " " << std::max(options_->spatialResolution.value(),
                                              0.1 * distanceSum) *
                            0.001 << "[km]" << std::endl;
    }
  }
  if (errorCategory == 0 || ((rejectedObservation->get().getObservationStatistics().
                              speedAveraged) >
                             options_->maxSpeed.value())) {
    oops::Log::trace() << "CheckShipTrack: cannot decide between station id " <<
                          trackId << " observations " <<
                          (observationAfterFastestSegment - 1)->get().getObservationNumber() <<
                          " " << observationAfterFastestSegment->get().getObservationNumber() <<
                          " rejecting both." << std::endl;
    errorCategory += 100;
    if (options_->testingMode.value() && firstIterativeRemoval) {
      std::vector<size_t> observationNumbersAroundFastest{
        (observationAfterFastestSegment - 1)->get().getObservationNumber(),
            observationAfterFastestSegment->get().getObservationNumber()};
      diagnostics_->storeFirstIterativeRemovalInfo(
            std::make_pair(observationNumbersAroundFastest, errorCategory));
    }
    if (rejectedObservation == observationAfterFastestSegment) {
      fail(observationAfterFastestSegment - 1);
    } else {
      fail(observationAfterFastestSegment);
    }
    track.erase(observationAfterFastestSegment - 1, observationAfterFastestSegment);
  } else {
    if (options_->testingMode.value() && firstIterativeRemoval) {
      std::vector<size_t> rejectedObservationNumber{rejectedObservation->get().
            getObservationNumber()};
      diagnostics_->storeFirstIterativeRemovalInfo(
            std::make_pair(rejectedObservationNumber,
                           errorCategory));
    }
    oops::Log::trace() << "CheckShipTrack: rejecting station " << trackId << " observation " <<
                          rejectedObservation->get().getObservationNumber() << "\n" <<
                          "Error category: " << errorCategory << "\n" <<
                          "rejection candidates: " <<
                          (observationAfterFastestSegment - 1)->get().getObservationNumber() <<
                          " " << observationAfterFastestSegment->get().getObservationNumber() <<
                          "\n" << "speeds: " << (observationAfterFastestSegment - 1)->get().
                          getObservationStatistics().speed << " " <<
                          observationAfterFastestSegment->get().getObservationStatistics().speed <<
                          "\n" << (observationAfterFastestSegment - 1)->get().
                          getObservationStatistics().angle << " " <<
                          observationAfterFastestSegment->get().
                          getObservationStatistics().angle << "\n";
    track.erase(rejectedObservation);
  }
}
/// \todo Trace output will need to be changed to match that of OPS (indices, LWin)

/// \brief Calculates all of the statistics that require only two
/// adjacent \p TrackObservations, storing within the righthand observation.
///
/// This includes
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
  if (firstIteration) {
    adjustTwoObservationStatistics(options);
  }
}

void TrackCheckShip::TrackObservation::resetObservationCalculations() {
  this->setDistance(0.0);
  this->setSpeed(0.0);
  this->setDistanceAveraged(0.0);
  this->setSpeedAveraged(0.0);
  this->setAngle(0.0);
}

/// Calculates all of the statistics that require three
/// consecutive \p TrackObservations,
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
/// \p calculateTwoObservationValues and \p calculateThreeObservationValues for the
/// non-edge-case observations.
void TrackCheckShip::calculateTrackSegmentProperties(
    const std::vector<std::reference_wrapper<TrackObservation>> &trackObservations,
    CalculationMethod calculationMethod) const {
  if (trackObservations.size()) {
    if (calculationMethod == MAINLOOP)
      trackObservations[0].get().resetObservationCalculations();
    for (size_t obsIdx = 1; obsIdx < trackObservations.size(); obsIdx++) {
      TrackObservation &obs = trackObservations[obsIdx].get();
      TrackObservation &prevObs = trackObservations[obsIdx - 1].get();
      obs.calculateTwoObservationValues(prevObs, calculationMethod == FIRSTITERATION, *options_);
      if (obsIdx > 1) {
        const TrackObservation &prevPrevObs = trackObservations[obsIdx - 2].get();
        prevObs.calculateThreeObservationValues(prevPrevObs, obs,
                                                calculationMethod == FIRSTITERATION, *options_);
      }
      if (calculationMethod == FIRSTITERATION && (obsIdx == trackObservations.size() - 1)) {
        int potentialDenominator = trackObservations.size() - 1 -
            obs.getFullTrackStatistics()->numShort_ - obs.getFullTrackStatistics()->numFast_;
        (obs.getFullTrackStatistics()->meanSpeed_) = (obs.getFullTrackStatistics()->sumSpeed_) /
            std::max(1, potentialDenominator);
      }
    }
    if (options_->testingMode.value() && calculationMethod != MAINLOOP) {
      std::vector<TrackCheckShip::ObservationStatistics> obsStats;
      for (size_t obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
        obsStats.push_back(trackObservations[obsIdx].get().getObservationStatistics());
      }
      auto trackStats = *(trackObservations[0].get().getFullTrackStatistics());
      if (calculationMethod == FIRSTITERATION)
        diagnostics_->storeInitialCalculationResults(std::make_pair(obsStats, trackStats));
    }
  }
}

/// \brief Keeps track of 0-distanced, short, and fast track segments,
/// as well as incrementing \p sumSpeed_ for normal track segments.
void TrackCheckShip::TrackObservation::adjustTwoObservationStatistics
(const TrackCheckShipParameters &options) const {
  util::Duration hour{"PT1H"};
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
