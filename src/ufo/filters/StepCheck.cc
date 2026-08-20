/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/StepCheck.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

StepCheck::StepCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                     ioda::ObsDataVector<int> & flags,
                     ioda::ObsDataVector<float> & obserr)
    : FilterBase(obsdb, parameters, flags, obserr), options_(parameters) {
  oops::Log::trace() << "StepCheck constructor" << std::endl;

  // Validate parameter combinations
  if (options_.percentageStepTolerance.value()) {
    if (options_.numberStepTolerance.value()) {
      throw eckit::UserError(R"(If percentage step tolerance is set,
            number step tolerance must not also be set.)",
                             Here());
    }
    if (options_.percentageStepTolerance.value().value() < 0 ||
        options_.percentageStepTolerance.value().value() > 100) {
      throw eckit::UserError(
          R"(Percentage step tolerance must be between 0 and 100.)", Here());
    }
    if (options_.useAverageStep.value()) {
      throw eckit::UserError(R"(Percentage step tolerance cannot be used with
            use average step.)",
                             Here());
    }
  }

  if (options_.useAverageStep.value()) {
    if (options_.numberStepTolerance.value() ||
        options_.percentageStepTolerance.value()) {
      throw eckit::UserError(R"(Use average step cannot be used with
            number step tolerance or percentage step tolerance.)",
                             Here());
    }
  }

  // Chunking requires use average step to be enabled
  if (options_.chunkSize.value() && !options_.useAverageStep.value()) {
    throw eckit::UserError(
        R"(Chunking cannot be used without use average step being set to true.)",
        Here());
  }

  oops::Log::debug() << "StepCheck: config = " << options_ << '\n';
}

StepCheck::~StepCheck() {
  oops::Log::trace() << "StepCheck destructor" << std::endl;
}

void StepCheck::applyFilter(const std::vector<bool> &apply,
                            const Variables &filtervars,
                            std::vector<std::vector<bool>> &flagged) const {
  oops::Log::trace() << "StepCheck applyFilter start" << std::endl;

  // Create ObsAccessor for handling grouped observations
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(
      options_.stationIdVariable, obsdb_, false);
  const std::vector<size_t> validObsIds =
      obsAccessor.getValidObservationIds(apply);

  // Split observations into independent groups (e.g., by station)
  RecursiveSplitter splitter =
      obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  TrackCheckUtils::sortTracksChronologically(validObsIds, obsAccessor,
                                             splitter);

  std::vector<bool> isRejected(obsAccessor.totalNumObservations(), false);
  const std::vector<std::string> filterVariables =
      filtervars.toOopsVariables().variables();

  // Process each variable
  for (const std::string &variable : filterVariables) {
    if (!obsdb_.has("ObsValue", variable)) {
      throw eckit::UserError(
          "StepCheck Error: ObsValue vector for " + variable + " not found.",
          Here());
    }

    const std::vector<float> variableValues =
        obsAccessor.getFloatVariableFromObsSpace("ObsValue", variable);

    // Process each station/group
    for (auto station : splitter.multiElementGroups()) {
      const std::vector<float> stationData = collectStationVariableData(
          station.begin(), station.end(), validObsIds, variableValues);

      processStation(station.begin(), station.end(), validObsIds, stationData,
                     isRejected);
    }
  }

  obsAccessor.flagRejectedObservations(isRejected, flagged);
  oops::Log::trace() << "StepCheck applyFilter complete" << std::endl;
}

void StepCheck::print(std::ostream &os) const {
  os << "StepCheck: config = " << options_ << '\n';
}

std::vector<float> StepCheck::collectStationVariableData(
    std::vector<size_t>::const_iterator stationObsIndicesBegin,
    std::vector<size_t>::const_iterator stationObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &globalData) const {
  std::vector<float> stationData;
  stationData.reserve(stationObsIndicesEnd - stationObsIndicesBegin);
  for (std::vector<size_t>::const_iterator it = stationObsIndicesBegin;
       it != stationObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds.at(*it);
    stationData.push_back(globalData[obsId]);
  }
  return stationData;
}

float StepCheck::calculateDifference(const float val1, const float val2) const {
  if (options_.circularPeriod.value()) {
    // Circular difference calculation
    const float period = options_.circularPeriod.value().value();
    float diff = std::abs(val2 - val1);
    // Take the shorter path around the circle
    if (diff > period / 2.0f) {
      diff = period - diff;
    }
    return diff;
  } else {
    // Standard absolute difference
    return std::abs(val2 - val1);
  }
}

void StepCheck::processStation(std::vector<size_t>::const_iterator stationBegin,
                               std::vector<size_t>::const_iterator stationEnd,
                               const std::vector<size_t> &validObsIds,
                               const std::vector<float> &stationData,
                               std::vector<bool> &isRejected) const {
  const float missingFloat = util::missingValue<float>();

  // Collect valid (non-missing) observations and their indices
  std::vector<float> validData;
  std::vector<size_t> validIndices;
  for (size_t i = 0; i < stationData.size(); ++i) {
    if (stationData[i] != missingFloat) {
      validData.push_back(stationData[i]);
      validIndices.push_back(i);
    }
  }

  if (validData.size() < 2) {
    return;  // Need at least 2 observations for step detection
  }

  // Determine which observations should be flagged
  const std::vector<size_t> indicesToFlag =
      determineObservationsToFlag(validData);

  // Flag the appropriate observations
  for (size_t idx : indicesToFlag) {
    const size_t stationLocalIdx = validIndices[idx];
    const size_t obsId = validObsIds.at(*(stationBegin + stationLocalIdx));
    isRejected[obsId] = true;
  }
}

std::vector<float> StepCheck::calculateConsecutiveSteps(
    const std::vector<float> &values) const {
  std::vector<float> steps;
  steps.reserve(values.size() - 1);
  for (size_t i = 1; i < values.size(); ++i) {
    steps.push_back(calculateDifference(values[i - 1], values[i]));
  }
  return steps;
}

bool StepCheck::exceedsThreshold(const float value,
                                 const float threshold) const {
  return options_.inclusiveStepThreshold.value() ? (value >= threshold)
                                                 : (value > threshold);
}

float StepCheck::calculateAverageStep(
    const std::vector<float> &validData) const {
  if (validData.size() < 2) {
    return -1.0f;  // Insufficient data
  }

  // No chunking: calculate mean of all consecutive steps
  if (!options_.chunkSize.value()) {
    return calculateMean(calculateConsecutiveSteps(validData));
  }

  // Chunking mode: split into chunks, calculate mean step per chunk,
  // then average across chunks
  const size_t chunkSize = options_.chunkSize.value().value();

  if (validData.size() < chunkSize) {
    return -1.0f;  // Not enough data for even one chunk
  }

  // Create chunks
  std::vector<std::vector<float>> chunks;
  for (size_t i = 0; i < validData.size(); i += chunkSize) {
    const size_t currentChunkSize = std::min(chunkSize, validData.size() - i);
    chunks.emplace_back(validData.begin() + i,
                        validData.begin() + i + currentChunkSize);
  }

  // Handle incomplete last chunk
  if (chunks.back().size() < chunkSize &&
      options_.ignoreLastChunkIfIncomplete.value()) {
    chunks.pop_back();
  }

  if (chunks.empty()) {
    return -1.0f;  // No chunks left
  }

  // Calculate mean step within each non-stuck chunk
  std::vector<float> chunkMeanSteps;
  for (const auto &chunk : chunks) {
    if (options_.removeStuckChunks.value() && isChunkStuck(chunk)) {
      continue;  // Skip stuck chunks
    }
    if (chunk.size() < 2) {
      continue;  // Need at least 2 values for step calculation
    }
    chunkMeanSteps.push_back(calculateMean(calculateConsecutiveSteps(chunk)));
  }

  if (chunkMeanSteps.empty()) {
    return -1.0f;  // All chunks were stuck or invalid
  }

  return calculateMean(chunkMeanSteps);
}

std::vector<size_t> StepCheck::determineObservationsToFlag(
    const std::vector<float> &values) const {
  if (values.size() < 2) {
    return {};
  }

  const float stepThreshold = options_.stepThreshold.value();

  // Average step mode (includes chunking case): flag all or nothing
  if (options_.useAverageStep.value()) {
    const float averageStep = calculateAverageStep(values);

    // Negative value indicates insufficient data (e.g., all chunks stuck)
    if (averageStep < 0.0f) {
      return {};
    }

    if (exceedsThreshold(averageStep, stepThreshold)) {
      std::vector<size_t> allIndices(values.size());
      std::iota(allIndices.begin(), allIndices.end(), 0);
      return allIndices;
    }
    return {};
  }

  // Individual step mode: find steps that exceed threshold
  std::vector<size_t> exceedingStepIndices;
  const std::vector<float> steps = calculateConsecutiveSteps(values);
  for (size_t i = 0; i < steps.size(); ++i) {
    if (exceedsThreshold(steps[i], stepThreshold)) {
      exceedingStepIndices.push_back(i + 1);  // Index of second value in pair
    }
  }

  // Apply tolerance logic to decide whether to flag
  if (!shouldFlagBasedOnTolerance(exceedingStepIndices.size(), steps.size())) {
    return {};
  }

  return exceedingStepIndices;
}

bool StepCheck::shouldFlagBasedOnTolerance(const size_t exceedingCount,
                                           const size_t totalSteps) const {
  if (options_.percentageStepTolerance.value()) {
    const float percentageTolerance =
        options_.percentageStepTolerance.value().value();
    const float exceedingPercentage = (exceedingCount * 100.0f) / totalSteps;
    return exceedingPercentage > percentageTolerance;
  }

  if (options_.numberStepTolerance.value()) {
    const size_t numberStepTolerance =
        options_.numberStepTolerance.value().value();
    return exceedingCount > numberStepTolerance;
  }

  // No tolerance: flag if any steps exceed
  return exceedingCount > 0;
}

bool StepCheck::isChunkStuck(const std::vector<float> &chunkValues) const {
  if (chunkValues.empty()) return true;

  const float tolerance = options_.chunkStuckTolerance.value();

  auto[minIt, maxIt] =
      std::minmax_element(chunkValues.begin(), chunkValues.end());
  float range = *maxIt - *minIt;
  if (options_.circularPeriod.value()) {
    const float period = options_.circularPeriod.value().value();
    range = std::min(range, period - range);
  }
  return range <= tolerance;
}

float StepCheck::calculateMean(const std::vector<float> &values) const {
  if (values.empty()) return 0.0f;
  const float sum = std::accumulate(values.begin(), values.end(), 0.0f);
  return sum / static_cast<float>(values.size());
}

}  // namespace ufo
