/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/StuckCheck.h"

#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/optional.hpp>
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
                        ioda::ObsDataVector<int> & flags,
                        ioda::ObsDataVector<float> & obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::trace() << "StuckCheck constructor" << std::endl;
  if (options_.core.percentageStuckTolerance.value()) {
    if ((options_.core.numberStuckTolerance.value()) ||
           (options_.core.timeStuckTolerance.value()) ||
           (options_.core.numberStuckToleranceVariable.value())) {
      throw eckit::UserError(R"(If percentage stuck tolerance is set,
        neither number stuck tolerance, number stuck tolerance variable,
        nor time stuck tolerance should be set.)", Here());
    }
    if ((options_.core.percentageStuckTolerance.value().value() < 0) ||
        (options_.core.percentageStuckTolerance.value().value() > 100)) {
      throw eckit::UserError(R"(Percentage stuck tolerance must be between 0 and 100.)", Here());
    }
    if (options_.core.minimumAllowedStuck.value() < 1) {
      throw eckit::UserError(R"(Minimum allowed stuck must be at least 1.)",
                             Here());
    }
  } else {  // no percentage stuck tolerance
    if (options_.core.numberStuckTolerance.value() &&
        options_.core.numberStuckToleranceVariable.value()) {
      throw eckit::UserError(R"(If number stuck tolerance variable is set,
            number stuck tolerance must not also be set.)",
                             Here());
    }
    if (!(options_.core.numberStuckTolerance.value() ||
          options_.core.numberStuckToleranceVariable.value()) ||
        !options_.core.timeStuckTolerance.value()) {
      throw eckit::UserError(R"(If percentage stuck tolerance is not set,
            then number stuck tolerance or number stuck tolerance variable,
            and time stuck tolerance must be set.)",
                             Here());
    }
  }
  oops::Log::debug() << "StuckCheck: config = " << options_ << '\n';
}

StuckCheck::~StuckCheck() {
  oops::Log::trace() << "StuckCheck destructor" << std::endl;
}

/// The filter removes observations if they are part of a 'streak'. A streak is where the number
/// of identical observation values in sequence (for a given variable) is greater than a user
/// defined count.
/// To be a streak, it must also continue for longer than a user-defined duration or every
/// observation in the station's group must have an identical value.
void StuckCheck::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "StuckCheck applyFilter start" << std::endl;
  // 3rd arg: recordsAreSingleObs = false for Stuck Check.
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(options_.stationIdVariable,
                                                               obsdb_,
                                                               false);
  const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply);
  obsGroupDateTimes_ = obsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "dateTime");
  if (options_.core.numberStuckToleranceVariable.value()) {
    const Variable &numStuckVar =
        *options_.core.numberStuckToleranceVariable.value();
    if (obsdb_.dtype(numStuckVar.group(), numStuckVar.variable()) !=
        ioda::ObsDtype::Integer) {
      throw eckit::UserError(
          "The number stuck tolerance variable must have integer type.",
          Here());
    }
    numberStuckToleranceVarValues_ = obsAccessor.getIntVariableFromObsSpace(
        numStuckVar.group(), numStuckVar.variable());
  }
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
      const boost::optional<size_t> stuckTolerance = getStuckToleranceForRecord(
          validObsIds, station.begin(), station.end());
      if (!stuckTolerance) {
        stationNumber++;
        continue;
      }
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
                                                  stationId,
                                                  *stuckTolerance);
            }
          } else {  // streak ended in the previous observation
            StuckCheck::potentiallyRejectStreak(station.begin(),
                                                station.end(),
                                                validObsIds,
                                                firstSameValueIndex,
                                                observationIndex - 1,
                                                isRejected,
                                                stationId,
                                                *stuckTolerance);
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
  oops::Log::trace() << "StuckCheck applyFilter complete" << std::endl;
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

/// Determines the number stuck tolerance to use for the given record,
/// considering all three tolerance modes (percentage, per-observation variable,
/// fixed). Returns boost::none if no rejection should be attempted for this
/// record (computed percentage tolerance is below minimumAllowedStuck, or all
/// tolerance variable values are missing).
boost::optional<size_t> StuckCheck::getStuckToleranceForRecord(
    const std::vector<size_t> &validObsIds,
    std::vector<size_t>::const_iterator stationObsIndicesBegin,
    std::vector<size_t>::const_iterator stationObsIndicesEnd) const {
  const size_t stationLength = stationObsIndicesEnd - stationObsIndicesBegin;
  const size_t startObsIndex = validObsIds.at(*stationObsIndicesBegin);
  if (options_.core.percentageStuckTolerance.value()) {
    const size_t percentageBaseCount =
        options_.core.percentageStuckToleranceBasedOnIntervals.value()
            ? stationLength - 1
            : stationLength;
    const size_t computed =
        std::round(options_.core.percentageStuckTolerance.value().value() *
                   percentageBaseCount / 100.0);
    if (computed <
        static_cast<size_t>(options_.core.minimumAllowedStuck.value()))
      return boost::none;
    return computed;
  }
  if (options_.core.numberStuckToleranceVariable.value()) {
    const int missingInt = util::missingValue<int>();
    boost::optional<int> foundTolerance;
    for (auto it = stationObsIndicesBegin; it != stationObsIndicesEnd; ++it) {
      const size_t obsIndex = validObsIds.at(*it);
      const int currentObsTolerance =
          numberStuckToleranceVarValues_.at(obsIndex);
      if (currentObsTolerance == missingInt) continue;
      if (!foundTolerance) {
        foundTolerance = currentObsTolerance;
      } else if (currentObsTolerance != *foundTolerance) {
        throw eckit::BadValue(
            "StuckCheck Error: Not all number stuck tolerance values are equal "
            "for record starting at index " +
                std::to_string(startObsIndex) + ".",
            Here());
      }
    }
    if (!foundTolerance) return boost::none;
    if (*foundTolerance < 0) {
      throw eckit::BadValue(
          "The number stuck tolerance variable must be non-negative. "
          "Invalid value found for record starting at index " +
              std::to_string(startObsIndex) + ".",
          Here());
    }
    return static_cast<size_t>(*foundTolerance);
  }
  return static_cast<size_t>(
      options_.core.numberStuckTolerance.value().value());
}

void StuckCheck::potentiallyRejectStreak(
    std::vector<size_t>::const_iterator stationIndicesBegin,
    std::vector<size_t>::const_iterator stationIndicesEnd,
    const std::vector<size_t> &validObsIds,
    size_t startOfStreakIndex,
    size_t endOfStreakIndex,
    std::vector<bool> &isRejected,
    std::string stationId,
    size_t stuckTolerance) const {

  auto getObservationTime = [this, &stationIndicesBegin, &validObsIds] (
      size_t offsetFromBeginning)->util::DateTime{
    const size_t obsIndex = validObsIds.at(*(stationIndicesBegin + offsetFromBeginning));
    return obsGroupDateTimes_.at(obsIndex);
  };

  auto rejectObservation = [&validObsIds, &isRejected, &stationIndicesBegin, &stationId](
      size_t observationIndex) {
    const size_t obsIndex = validObsIds.at(*(stationIndicesBegin + observationIndex));
    isRejected[obsIndex] = true;
  };

  const size_t streakLength = endOfStreakIndex - startOfStreakIndex;
  if (streakLength < stuckTolerance) {
    return;
  }

  if (options_.core.timeStuckTolerance.value()) {
    const size_t stationLength = stationIndicesEnd - stationIndicesBegin;
    if (streakLength < stationLength - 1) {
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
