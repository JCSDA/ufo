/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/StuckCheck.h"

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/none.hpp>
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/StuckCheckParameters.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

StuckCheck::StuckCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                        std::shared_ptr<ioda::ObsDataVector<int> > flags,
                        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  if (options_.core.percentageStuckTolerance.value()) {
    if ((options_.core.numberStuckTolerance.value()) ||
           (options_.core.timeStuckTolerance.value())) {  // must NOT be set if percentage set
      throw eckit::UserError(R"(If percentage stuck tolerance is set,
            neither number stuck tolerance nor time stuck tolerance should be set.)", Here());
    }
    if ((options_.core.percentageStuckTolerance.value().value() < 0) ||
        (options_.core.percentageStuckTolerance.value().value() > 100)) {  // must be 0-100
      throw eckit::UserError(R"(Percentage stuck tolerance must be between 0 and 100.)", Here());
    }
  } else {  // no percentage stuck tolerance, so others MUST be set
    if (!options_.core.numberStuckTolerance.value() || !options_.core.timeStuckTolerance.value()) {
      throw eckit::UserError(R"(If percentage stuck tolerance is not set,
            then both number stuck tolerance and time stuck tolerance must be set.)", Here());
    }
  }
  obsGroupDateTimes_.reset(new std::vector<util::DateTime>);
  oops::Log::debug() << "StuckCheck: config = " << options_ << '\n';
}

// Required for the correct destruction of ObsGroupDateTimes_.
StuckCheck::~StuckCheck()
{}

/// The filter removes observations if they are part of a 'streak'. A streak is where the number
/// of identical observation values in sequence (for a given variable) is greater than a user
/// defined count.
/// To be a streak, it must also continue for longer than a user-defined duration or every
/// observation in the station's group must have an identical value.
void StuckCheck::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  // 3rd arg: recordsAreSingleObs = false for Stuck Check.
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(options_.stationIdVariable,
                                                               obsdb_,
                                                               false);
  const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply);
  *obsGroupDateTimes_ = obsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "dateTime");
  // Create groups based on record number (assumed station ID) or category variable
  // (stationIdVariable) or otherwise assume observations all taken by the same station (1 group)
  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  TrackCheckUtils::sortTracksChronologically(validObsIds, obsAccessor, splitter);
  std::vector<bool> isRejected(obsAccessor.totalNumObservations(), false);
  std::vector<std::string> filterVariables = filtervars.toOopsVariables().variables();
  // Iterates through observations to see how long each variable is stuck on one observation
  for (std::string const& variable : filterVariables) {
    size_t stationNumber = 0;
    if (!obsdb_.has("ObsValue", variable)) {
      std::string errorMessage =
          "StuckCheck Error: ObsValue vector for " + variable + " not found.\n";
      throw std::invalid_argument(errorMessage);
    }
    const std::vector<float> variableValues = obsAccessor.getFloatVariableFromObsSpace(
          "ObsValue", variable);
    const float missingFloat = util::missingValue<float>();
    for (auto station : splitter.multiElementGroups()) {
      std::string stationId = std::to_string(stationNumber);
      std::vector<float> variableDataStation = collectStationVariableData(
            station.begin(), station.end(), validObsIds, variableValues);
      // the working variable's value associated with the prior observation
      float previousObservationValue;
      float currentObservationValue;
      size_t firstSameValueIndex = 0;  // the first observation in the current streak
      for (size_t observationIndex = 0; observationIndex < variableDataStation.size();
           observationIndex++) {
        currentObservationValue = variableDataStation.at(observationIndex);
        if (currentObservationValue == missingFloat) {
          continue;
        }
        if (observationIndex == 0) {
          previousObservationValue = currentObservationValue;
        } else {
          if (currentObservationValue == previousObservationValue) {
            // If the last observation of the track is part of a streak, the full streak will need
            // to be checked at this point.
            if (observationIndex == variableDataStation.size() - 1) {
              StuckCheck::potentiallyRejectStreak(station.begin(),
                                                  station.end(),
                                                  validObsIds,
                                                  firstSameValueIndex,
                                                  observationIndex,
                                                  isRejected,
                                                  stationId);
            }
          } else {  // streak ended in the previous observation
            StuckCheck::potentiallyRejectStreak(station.begin(),
                                                station.end(),
                                                validObsIds,
                                                firstSameValueIndex,
                                                observationIndex - 1,
                                                isRejected,
                                                stationId);
            // start the streak with the current observation and reset the count to 1
            firstSameValueIndex = observationIndex;
            previousObservationValue = currentObservationValue;
          }
        }
      }
      stationNumber++;
    }
  }
  obsAccessor.flagRejectedObservations(isRejected, flagged);
}

void StuckCheck::print(std::ostream & os) const {
  os << "StuckCheck: config = " << options_ << '\n';
}

/// \returns a vector containing all of the necessary data to run this filter for each observation,
/// stored by observation.
std::vector<float> StuckCheck::collectStationVariableData(
    std::vector<size_t>::const_iterator stationObsIndicesBegin,
    std::vector<size_t>::const_iterator stationObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &globalData) const {
  std::vector<float> stationData;
  stationData.reserve(stationObsIndicesEnd - stationObsIndicesBegin);
  size_t observationNumber = 0;
  for (std::vector<size_t>::const_iterator it = stationObsIndicesBegin;
       it != stationObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds.at(*it);
    stationData.push_back(globalData[obsId]);
    observationNumber++;
  }
  return stationData;
}

void StuckCheck::potentiallyRejectStreak(
    std::vector<size_t>::const_iterator stationIndicesBegin,
    std::vector<size_t>::const_iterator stationIndicesEnd,
    const std::vector<size_t> &validObsIds,
    size_t startOfStreakIndex,
    size_t endOfStreakIndex,
    std::vector<bool> &isRejected,
    std::string stationId = "") const {

  auto getObservationTime = [this, &stationIndicesBegin, &validObsIds] (
      size_t offsetFromBeginning)->util::DateTime{
    const size_t obsIndex = validObsIds.at(*(stationIndicesBegin + offsetFromBeginning));
    return obsGroupDateTimes_->at(obsIndex);
  };

  auto rejectObservation = [&validObsIds, &isRejected, &stationIndicesBegin, &stationId](
      size_t observationIndex) {
    const size_t obsIndex = validObsIds.at(*(stationIndicesBegin + observationIndex));
    isRejected[obsIndex] = true;
  };

  const size_t streakLength = endOfStreakIndex - startOfStreakIndex + 1;
  const size_t stationLength = stationIndicesEnd - stationIndicesBegin;
  size_t numberStuckTolerance;
  if (options_.core.percentageStuckTolerance.value()) {
    numberStuckTolerance =
        std::round(options_.core.percentageStuckTolerance.value().value()*stationLength/100.0);
    if (numberStuckTolerance < 2) {
      return;
    }
  } else {
    numberStuckTolerance = options_.core.numberStuckTolerance.value().value();
  }
  if (streakLength <= numberStuckTolerance) {
      return;
  }

  if (!(options_.core.percentageStuckTolerance.value())) {
    if (streakLength < stationLength) {
      const util::DateTime firstStreakObservationTime = getObservationTime(startOfStreakIndex);
      const util::DateTime lastStreakObservationTime = getObservationTime(endOfStreakIndex);
      const util::Duration streakDuration = lastStreakObservationTime - firstStreakObservationTime;
      if (streakDuration <= options_.core.timeStuckTolerance.value().value()) {
        return;
      }
    }
  }
  // reject all observations within streak
  for (size_t indexToReject = startOfStreakIndex;
       indexToReject <= endOfStreakIndex;
       indexToReject++)
    rejectObservation(indexToReject);
}

}  // namespace ufo
