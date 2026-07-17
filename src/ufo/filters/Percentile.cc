/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/Percentile.h"

#include <algorithm>
#include <cmath>
#include <ostream>
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

Percentile::Percentile(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                       ioda::ObsDataVector<int> &flags,
                       ioda::ObsDataVector<float> &obserr)
    : FilterBase(obsdb, parameters, flags, obserr), options_(parameters) {
  oops::Log::trace() << "Percentile constructor" << std::endl;

  // Validate parameters
  if (!options_.lowerPercentiles.value() &&
      !options_.upperPercentiles.value()) {
    throw eckit::BadParameter(
        "At least one of 'lower percentiles' or 'upper percentiles' must be "
        "specified",
        Here());
  }
  oops::Log::debug() << "Percentile: config = " << options_ << '\n';
}

Percentile::~Percentile() {
  oops::Log::trace() << "Percentile destructor" << std::endl;
}

void Percentile::applyFilter(const std::vector<bool> &apply,
                             const Variables &filtervars,
                             std::vector<std::vector<bool>> &flagged) const {
  oops::Log::trace() << "Percentile applyFilter start" << std::endl;

  // Parse parameters - expand single values to match number of variables
  const size_t nlocs = obsdb_.nlocs();
  const std::vector<std::string> filterVariables =
      filtervars.toOopsVariables().variables();
  const size_t nvars = filterVariables.size();
  std::vector<float> lowerPercentiles;
  std::vector<float> upperPercentiles;
  std::vector<bool> inclusiveCentralRanges =
      options_.inclusiveCentralRange.value();

  // Validate lower percentiles
  if (options_.lowerPercentiles.value()) {
    lowerPercentiles = *options_.lowerPercentiles.value();
    if (lowerPercentiles.size() != nvars) {
      throw eckit::BadParameter(
          "Number of lower percentiles must match number of filter variables",
          Here());
    }
  } else {
    lowerPercentiles.resize(nvars, 0.0f);
  }

  // Validate upper percentiles
  if (options_.upperPercentiles.value()) {
    upperPercentiles = *options_.upperPercentiles.value();
    if (upperPercentiles.size() != nvars) {
      throw eckit::BadParameter(
          "Number of upper percentiles must match number of filter variables",
          Here());
    }
  } else {
    upperPercentiles.resize(nvars, 100.0f);
  }

  // Expand default single inclusive value to all variables
  if (inclusiveCentralRanges.size() == 1 && nvars > 1) {
    inclusiveCentralRanges.resize(nvars, inclusiveCentralRanges.front());
  }
  // Validate inclusive central range
  if (inclusiveCentralRanges.size() != nvars) {
    throw eckit::BadParameter(
        "Number of inclusive central range values must match number of filter "
        "variables",
        Here());
  }

  // Validate percentile ranges
  for (size_t i = 0; i < nvars; ++i) {
    if (lowerPercentiles[i] < 0.0f || lowerPercentiles[i] > 100.0f ||
        upperPercentiles[i] < 0.0f || upperPercentiles[i] > 100.0f) {
      throw eckit::BadParameter("Percentile values must be in range [0, 100]",
                                Here());
    }
    if (lowerPercentiles[i] > upperPercentiles[i]) {
      throw eckit::BadParameter("Lower percentile must be <= upper percentiles",
                                Here());
    }
  }

  const bool stationIdVariableUsed =
      options_.stationIdVariable.value() != boost::none;
  const bool recordsAreGrouped = !obsdb_.obs_group_vars().empty();
  if (stationIdVariableUsed && recordsAreGrouped) {
    throw eckit::BadParameter(
        "Using station_id_variable with record grouping is not supported "
        "for Percentile filter",
        Here());
  }
  if (!recordsAreGrouped && !stationIdVariableUsed) {
    throw eckit::BadParameter(
        "station_id_variable must be used when record grouping is not set "
        "for Percentile filter",
        Here());
  }

  // Create ObsAccessor for handling grouped observations - if an obs grouping
  // is defined in ObsSpace, use it; otherwise, use station_id_variable if
  // defined in the filter parameters.
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(
      options_.stationIdVariable, obsdb_, false);
  const std::vector<size_t> validObsIds =
      obsAccessor.getValidObservationIds(apply);

  // Split observations into independent groups (e.g., by station)
  RecursiveSplitter splitter =
      obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);

  // Vector to track rejected observations which can be indexed with
  // global observation indices (those returned by getValidObservationIds).
  std::vector<bool> isRejected(obsAccessor.totalNumObservations(), false);

  // Process each variable
  for (size_t varIndex = 0; varIndex < filterVariables.size(); ++varIndex) {
    const std::string &variable = filterVariables[varIndex];
    oops::Log::debug() << "Processing variable: " << variable << std::endl;

    // Get variable data
    const std::vector<float> variableData =
        obsAccessor.getFloatVariableFromObsSpace("ObsValue", variable);

    // Prepare output data vector - size to match variableData (global-sized
    // when using station_id_variable)
    std::vector<float> dataOutput(variableData.size(),
                                  util::missingValue<float>());

    const float lowerPercentile = lowerPercentiles[varIndex];
    const float upperPercentile = upperPercentiles[varIndex];
    const bool inclusive = inclusiveCentralRanges[varIndex];

    // Process each record (station/group)
    for (const auto &record : splitter.groups()) {
      processRecord(record, validObsIds, variableData,
                    lowerPercentile, upperPercentile, inclusive, isRejected,
                    dataOutput);
    }

    // Extract local portion of dataOutput for writing to obsdb
    // When observations are grouped into records, dataOutput is already
    // local-sized and correctly indexed.
    // When instead using station_id_variable without record grouping,
    // dataOutput is global-sized and needs distribution mapping
    std::vector<float> localDataOutput(nlocs, util::missingValue<float>());

    if (recordsAreGrouped) {
      // Observations are grouped into records - data is already local
      localDataOutput = dataOutput;
    } else {
      // Using station_id_variable without record grouping - need to extract
      // local observations from global data
      for (size_t localObsId = 0; localObsId < static_cast<size_t>(nlocs);
           ++localObsId) {
        const size_t globalObsId =
            obsdb_.distribution()->globalUniqueConsecutiveLocationIndex(
                localObsId);
        if (globalObsId < dataOutput.size()) {
          localDataOutput[localObsId] = dataOutput[globalObsId];
        }
      }
    }

    // Write output
    obsdb_.put_db("DerivedObsValue", variable, localDataOutput);
  }

  obsAccessor.flagRejectedObservations(isRejected, flagged);
  oops::Log::trace() << "Percentile applyFilter complete" << std::endl;
}

void Percentile::print(std::ostream &os) const {
  os << "Percentile: config = " << options_ << '\n';
}

void Percentile::processRecord(
    const RecursiveSplitter::Group &record,
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &variableData,
    const float lowerPercentile, const float upperPercentile,
    const bool inclusive, std::vector<bool> &isRejected,
    std::vector<float> &outputData) const {
  const float missing = util::missingValue<float>();
  const auto recordBegin = record.begin();
  const auto recordEnd = record.end();

  // Collect data for percentile calculation, making sure to not use missing values
  std::vector<float> recordData;
  recordData.reserve(recordEnd - recordBegin);
  for (auto it = recordBegin; it != recordEnd; ++it) {
    const size_t obsId = validObsIds.at(*it);
    const float value = variableData[obsId];
    if (value != missing) {
      recordData.push_back(value);
    }
  }

  if (recordData.empty()) {
    // No valid data to filter in this record
    return;
  }

  // Calculate percentile threshold values using linear interpolation between
  // closest datapoints where necessary to match the requested percentile values
  std::vector<float> sortedData = recordData;
  std::sort(sortedData.begin(), sortedData.end());
  const float lowerThreshold = calculatePercentileValue(sortedData, lowerPercentile);
  const float upperThreshold = calculatePercentileValue(sortedData, upperPercentile);

  for (auto it = recordBegin; it != recordEnd; ++it) {
    const size_t obsId = validObsIds.at(*it);
    const float value = variableData[obsId];

    // Skip missing values - still need to set outputData to missing but don't mark as rejected
    if (value == missing) {
      outputData[obsId] = missing;
      continue;
    }

    // Check if value is in the central percentile range
    const bool passes =
        inclusive ? (value >= lowerThreshold && value <= upperThreshold)
                  : (value > lowerThreshold && value < upperThreshold);

    // Update isRejected using global obsId
    if (!passes) {
      isRejected[obsId] = true;
    }

    // Write to output data using global obsId
    outputData[obsId] = passes ? value : missing;
  }
}

float Percentile::calculatePercentileValue(const std::vector<float> &sortedData,
                                           const float percentile) const {
  if (sortedData.empty()) {
    return util::missingValue<float>();
  }
  const float impliedIndex =
      (percentile / 100.0f) * static_cast<float>(sortedData.size() - 1);
  const size_t lowerIndex = static_cast<size_t>(std::floor(impliedIndex));
  const size_t upperIndex = static_cast<size_t>(std::ceil(impliedIndex));

  // If impliedIndex is an integer (lowerIndex == upperIndex) or upperIndex
  // is out of bounds, return the value at lowerIndex
  if (lowerIndex == upperIndex || upperIndex >= sortedData.size()) {
    return sortedData[lowerIndex];
  }
  // Otherwise linearly interpolate between the values at lowerIndex and
  // upperIndex
  const float indexFraction = impliedIndex - static_cast<float>(lowerIndex);
  return sortedData[lowerIndex] +
         indexFraction * (sortedData[upperIndex] - sortedData[lowerIndex]);
}

}  // namespace ufo
