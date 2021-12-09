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
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(options_.stationIdVariable, obsdb_);
  const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply);
  *obsGroupDateTimes_ = obsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "dateTime");
  // Create groups based on record number (assumed station ID) or category variable
  // (stationIdVariable) or otherwise assume observations all taken by the same station (1 group)
  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  TrackCheckUtils::sortTracksChronologically(validObsIds, obsAccessor, splitter);
  std::vector<bool> isRejected(validObsIds.size(), false);
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

  size_t streakLength = endOfStreakIndex - startOfStreakIndex + 1;
  if (streakLength <= options_.core.numberStuckTolerance) {
    return;
  }

  size_t stationLength = stationIndicesEnd - stationIndicesBegin;

  if (streakLength < stationLength) {
    util::DateTime firstStreakObservationTime = getObservationTime(startOfStreakIndex);
    util::DateTime lastStreakObservationTime = getObservationTime(endOfStreakIndex);
    util::Duration streakDuration = lastStreakObservationTime - firstStreakObservationTime;
    if (streakDuration <= options_.core.timeStuckTolerance) {
      return;
    }
  }
  // reject all observations within streak
  for (size_t indexToReject = startOfStreakIndex;
       indexToReject <= endOfStreakIndex;
       indexToReject++)
    rejectObservation(indexToReject);
}

}  // namespace ufo
