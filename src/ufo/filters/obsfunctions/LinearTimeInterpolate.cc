/*
 * (C) British Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/LinearTimeInterpolate.h"

#include <algorithm>
#include <regex>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<LinearTimeInterpolate> makerLinearTimeInterpolate_(
    "LinearTimeInterpolate");

// -----------------------------------------------------------------------------

LinearTimeInterpolate::LinearTimeInterpolate(
    const eckit::LocalConfiguration& conf)
    : invars_(),
      timestampsAreLiteral_(false),
      literalTimestamps_(),
      timestampsAreTimeOnly_(false) {
  oops::Log::trace() << "LinearTimeInterpolate constructor" << std::endl;

  options_.validateAndDeserialize(conf);

  const size_t numVars = options_.inputVariables.value().size();
  const auto& inputTimestamps = options_.inputTimestamps.value();
  const size_t numTimes = inputTimestamps.size();

  if (numVars < 2) {
    throw eckit::BadParameter(
        "LinearTimeInterpolate requires at least 2 input variables", Here());
  }
  if (numTimes < 2) {
    throw eckit::BadParameter(
        "LinearTimeInterpolate requires at least 2 input timestamps", Here());
  }
  if (numVars != numTimes) {
    throw eckit::BadParameter(
        "LinearTimeInterpolate requires the same number of input variables "
        "and timestamps",
        Here());
  }

  invars_ += options_.targetDateTime.value();
  for (const Variable& var : options_.inputVariables.value()) {
    invars_ += var;
  }

  const auto isVariableReference = [](const std::string& token) {
    return token.find('/') != std::string::npos;
  };

  const auto throwInvalidTimestampToken = [](const std::string& token) {
    throw eckit::BadParameter(
        "LinearTimeInterpolate: input timestamp '" + token +
            "' is neither a supported literal timestamp nor a variable "
            "reference of the form Group/variable",
        Here());
  };

  // Determine mode from the first entry, then enforce all others match
  util::DateTime firstParsed(1970, 1, 1, 0, 0, 0);
  bool firstIsTimeOnly = false;
  const bool firstIsLiteral = tryParseLiteralTimestamp(
      inputTimestamps[0], firstParsed, firstIsTimeOnly);

  if (!firstIsLiteral && !isVariableReference(inputTimestamps[0])) {
    throwInvalidTimestampToken(inputTimestamps[0]);
  }

  timestampsAreLiteral_ = firstIsLiteral;
  timestampsAreTimeOnly_ = firstIsTimeOnly;

  if (timestampsAreLiteral_) {
    literalTimestamps_.resize(numTimes);
    literalTimestamps_[0] = firstParsed;
    for (size_t i = 1; i < numTimes; ++i) {
      util::DateTime parsedLiteral;
      bool parsedLiteralIsTimeOnly = false;
      if (!tryParseLiteralTimestamp(inputTimestamps[i], parsedLiteral,
                                    parsedLiteralIsTimeOnly)) {
        throw eckit::BadParameter(
            "LinearTimeInterpolate: cannot mix literal timestamps with "
            "variable references in input timestamps",
            Here());
      }
      if (parsedLiteralIsTimeOnly != timestampsAreTimeOnly_) {
        throw eckit::BadParameter(
            "LinearTimeInterpolate: all input timestamps must use the same "
            "format (either all time-only THH:MM:SSZ or all full datetime "
            "YYYY-MM-DDTHH:MM:SSZ)",
            Here());
      }
      literalTimestamps_[i] = parsedLiteral;
    }
  } else {
    invars_ += Variable(inputTimestamps[0]);
    for (size_t i = 1; i < numTimes; ++i) {
      util::DateTime parsedLiteralCandidate;
      bool parsedIsTimeOnly = false;
      if (tryParseLiteralTimestamp(inputTimestamps[i], parsedLiteralCandidate,
                                   parsedIsTimeOnly)) {
        throw eckit::BadParameter(
            "LinearTimeInterpolate: cannot mix variable references with "
            "literal timestamps in input timestamps",
            Here());
      }
      if (!isVariableReference(inputTimestamps[i])) {
        throwInvalidTimestampToken(inputTimestamps[i]);
      }
      invars_ += Variable(inputTimestamps[i]);
    }
  }
}

// -----------------------------------------------------------------------------

LinearTimeInterpolate::~LinearTimeInterpolate() {
  oops::Log::trace() << "LinearTimeInterpolate destructor" << std::endl;
}

// -----------------------------------------------------------------------------

bool LinearTimeInterpolate::tryParseLiteralTimestamp(
    const std::string& varToken, util::DateTime& parsedDateTime,
    bool& parsedIsTimeOnly) const {
  static const std::regex timeOnlyPattern(R"(^T(\d{2}):(\d{2}):(\d{2})Z$)");
  static const std::regex fullDateTimePattern(
      R"(^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z$)");

  std::smatch match;
  if (std::regex_match(varToken, match, timeOnlyPattern)) {
    const int hour = std::stoi(match[1].str());
    const int minute = std::stoi(match[2].str());
    const int second = std::stoi(match[3].str());

    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 ||
        second > 59) {
      throw eckit::BadParameter(
          "LinearTimeInterpolate: Invalid time-only timestamp '" + varToken +
              "'. Expected THH:MM:SSZ with valid hour/minute/second values",
          Here());
    }

    // Use epoch date (1970-01-01) for time-only literals; date will be ignored
    parsedDateTime = util::DateTime(1970, 1, 1, hour, minute, second);
    parsedIsTimeOnly = true;
    return true;
  }

  if (std::regex_match(varToken, match, fullDateTimePattern)) {
    const int year = std::stoi(match[1].str());
    const int month = std::stoi(match[2].str());
    const int day = std::stoi(match[3].str());
    const int hour = std::stoi(match[4].str());
    const int minute = std::stoi(match[5].str());
    const int second = std::stoi(match[6].str());

    if (month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 ||
        hour > 23 || minute < 0 || minute > 59 || second < 0 || second > 59) {
      throw eckit::BadParameter(
          "LinearTimeInterpolate: Invalid full datetime timestamp '" +
              varToken +
              "'. Expected YYYY-MM-DDTHH:MM:SSZ with valid month/day/hour/"
              "minute/second values",
          Here());
    }

    parsedDateTime = util::DateTime(year, month, day, hour, minute, second);
    parsedIsTimeOnly = false;
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------

void LinearTimeInterpolate::compute(const ObsFilterData& in,
                                    ioda::ObsDataVector<float>& out) const {
  oops::Log::trace() << "LinearTimeInterpolate compute start" << std::endl;

  const size_t nlocs = in.nlocs();
  const size_t numPoints = options_.inputVariables.value().size();

  std::vector<util::DateTime> targetDateTimes(nlocs);
  in.get(options_.targetDateTime.value(), targetDateTimes);

  std::vector<ioda::ObsDataVector<float>> inputVariablesValues;
  inputVariablesValues.reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    const Variable& var = options_.inputVariables.value()[i];
    inputVariablesValues.emplace_back(in.obsspace(), var.toOopsObsVariables());
    in.get(var, inputVariablesValues[i]);
  }
  ASSERT(inputVariablesValues[0].nvars() == out.nvars());

  std::vector<ioda::ObsDataVector<util::DateTime>> inputDateTimes;
  if (!timestampsAreLiteral_) {
    inputDateTimes.reserve(numPoints);
    for (const std::string& timestampVariableName :
         options_.inputTimestamps.value()) {
      const Variable timestampVar(timestampVariableName);
      inputDateTimes.emplace_back(in.obsspace(),
                                  timestampVar.toOopsObsVariables());
      in.get(timestampVar, inputDateTimes.back());
    }
  }

  performInterpolation(out, targetDateTimes, nlocs, numPoints,
                       inputVariablesValues, inputDateTimes);

  oops::Log::trace() << "LinearTimeInterpolate compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

int64_t LinearTimeInterpolate::calculateTimeOffset(
    const util::DateTime& targetDateTime, const util::DateTime& inputDateTime,
    bool isTimeOnly) const {
  if (isTimeOnly) {
    // Ignore date components when calculating time offset for time-only format
    int targetHour, targetMinute, targetSecond;
    int inputHour, inputMinute, inputSecond;
    int dummyYear, dummyMonth, dummyDay;

    targetDateTime.toYYYYMMDDhhmmss(dummyYear, dummyMonth, dummyDay, targetHour,
                                    targetMinute, targetSecond);
    inputDateTime.toYYYYMMDDhhmmss(dummyYear, dummyMonth, dummyDay, inputHour,
                                   inputMinute, inputSecond);

    const int targetSeconds =
        targetHour * 3600 + targetMinute * 60 + targetSecond;
    const int inputSeconds = inputHour * 3600 + inputMinute * 60 + inputSecond;

    // Wrap to (-12h, +12h] for shortest path around the 24-hour clock from
    // target to input
    int64_t diff = inputSeconds - targetSeconds;
    constexpr int64_t halfDay = 43200;  // 12 hours in seconds
    constexpr int64_t fullDay = 86400;  // 24 hours in seconds
    if (diff > halfDay) {
      diff -= fullDay;
    } else if (diff <= -halfDay) {
      diff += fullDay;
    }
    return diff;
  } else {
    const util::Duration dt = inputDateTime - targetDateTime;
    return dt.toSeconds();
  }
}

// -----------------------------------------------------------------------------

std::vector<LinearTimeInterpolate::InterpolationPoint>
LinearTimeInterpolate::buildInterpolationPoints(
    size_t numPoints, size_t ichan, size_t iloc,
    const util::DateTime& targetDateTime,
    const std::vector<ioda::ObsDataVector<float>>& inputVariablesValues,
    const std::vector<ioda::ObsDataVector<util::DateTime>>& inputDateTimes)
    const {
  std::vector<InterpolationPoint> candidatePoints;
  candidatePoints.reserve(numPoints);

  for (size_t inputVarIndex = 0; inputVarIndex < numPoints; ++inputVarIndex) {
    util::DateTime inputDateTime(1970, 1, 1, 0, 0, 0);
    if (timestampsAreLiteral_) {
      inputDateTime = literalTimestamps_[inputVarIndex];
    } else {
      inputDateTime = inputDateTimes[inputVarIndex][0][iloc];
    }

    const int64_t timeOffset = calculateTimeOffset(
        targetDateTime, inputDateTime, timestampsAreTimeOnly_);

    candidatePoints.emplace_back(
        timeOffset, inputVariablesValues[inputVarIndex][ichan][iloc],
        inputVarIndex);
  }

  return candidatePoints;
}

// -----------------------------------------------------------------------------

void LinearTimeInterpolate::performInterpolation(
    ioda::ObsDataVector<float>& out,
    const std::vector<util::DateTime>& targetDateTimes, size_t nlocs,
    size_t numPoints,
    const std::vector<ioda::ObsDataVector<float>>& inputVariablesValues,
    const std::vector<ioda::ObsDataVector<util::DateTime>>& inputDateTimes)
    const {
  const float missing = util::missingValue<float>();

  for (size_t ichan = 0; ichan < out.nvars(); ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      const util::DateTime& targetDateTime = targetDateTimes[iloc];

      std::vector<InterpolationPoint> candidatePoints =
          buildInterpolationPoints(numPoints, ichan, iloc, targetDateTime,
                                   inputVariablesValues, inputDateTimes);

      std::vector<InterpolationPoint> validPoints =
          extractValidPoints(candidatePoints);

      if (validPoints.size() < 2) {
        out[ichan][iloc] = missing;
        continue;
      }

      std::sort(validPoints.begin(), validPoints.end(),
                [](const InterpolationPoint& a, const InterpolationPoint& b) {
                  return a.timeOffset < b.timeOffset;
                });

      std::vector<InterpolationPoint> uniquePoints =
          deduplicatePoints(validPoints, iloc, ichan);

      if (uniquePoints.empty()) {
        out[ichan][iloc] = missing;
        continue;
      }

      if (uniquePoints.size() < 2) {
        out[ichan][iloc] = missing;
        continue;
      }

      size_t initialUniquePointsIdx = 0;
      size_t subsequentUniquePointsIdx = 1;
      selectBracketingOrNearestPoints(uniquePoints, initialUniquePointsIdx,
                                subsequentUniquePointsIdx);

      if (!options_.allowGapInterpolation.value() &&
          hasGapInRange(candidatePoints, uniquePoints, initialUniquePointsIdx,
                        subsequentUniquePointsIdx)) {
        out[ichan][iloc] = missing;
        continue;
      }

      out[ichan][iloc] = evaluateLinearFitAtTargetTime(
          uniquePoints[initialUniquePointsIdx],
          uniquePoints[subsequentUniquePointsIdx]);
    }
  }
}

// -----------------------------------------------------------------------------

std::vector<LinearTimeInterpolate::InterpolationPoint>
LinearTimeInterpolate::extractValidPoints(
    const std::vector<InterpolationPoint>& candidatePoints) const {
  std::vector<InterpolationPoint> validPoints;
  validPoints.reserve(candidatePoints.size());

  for (const auto& point : candidatePoints) {
    if (point.isValid()) {
      validPoints.push_back(point);
    }
  }

  return validPoints;
}

// -----------------------------------------------------------------------------

std::vector<LinearTimeInterpolate::InterpolationPoint>
LinearTimeInterpolate::deduplicatePoints(
    const std::vector<InterpolationPoint>& validPoints, size_t iloc,
    size_t ichan) const {
  std::vector<InterpolationPoint> uniquePoints;
  uniquePoints.reserve(validPoints.size());
  uniquePoints.push_back(validPoints[0]);

  for (size_t i = 1; i < validPoints.size(); ++i) {
    if (validPoints[i].timeOffset == validPoints[i - 1].timeOffset) {
      // Duplicate timestamp found - check if values are also identical
      if (validPoints[i].value != validPoints[i - 1].value) {
        oops::Log::warning()
            << "LinearTimeInterpolate: duplicate timestamps with different "
               "values "
            << "detected at location " << iloc << ", channel " << ichan
            << " (indices " << validPoints[i - 1].inputVarIndex << " and "
            << validPoints[i].inputVarIndex << ")." << std::endl;
        return {};  // Return empty vector to signal conflict
      }
      // Same timestamp, same value - skip duplicate
    } else {
      // Different timestamp - add to unique points
      uniquePoints.push_back(validPoints[i]);
    }
  }

  return uniquePoints;
}

// -----------------------------------------------------------------------------

void LinearTimeInterpolate::selectBracketingOrNearestPoints(
    const std::vector<InterpolationPoint>& uniquePoints, size_t& idx1,
    size_t& idx2) const {
  if (uniquePoints[0].timeOffset >= 0) {
    // Target time is before the first point, so use the first two points for
    // extrapolation
    idx1 = 0;
    idx2 = 1;
  } else if (uniquePoints[uniquePoints.size() - 1].timeOffset <= 0) {
    // Target time is after the last point, so use the last two points for
    // extrapolation
    idx1 = uniquePoints.size() - 2;
    idx2 = uniquePoints.size() - 1;
  } else {
    // Target time is between points, find the two points that bracket the
    // target time for interpolation
    for (size_t i = 0; i < uniquePoints.size() - 1; ++i) {
      if (uniquePoints[i].timeOffset <= 0 &&
          uniquePoints[i + 1].timeOffset >= 0) {
        idx1 = i;
        idx2 = i + 1;
        break;
      }
    }
  }
}

// -----------------------------------------------------------------------------

float LinearTimeInterpolate::evaluateLinearFitAtTargetTime(
    const InterpolationPoint& point1, const InterpolationPoint& point2) const {
  const double t1 = static_cast<double>(point1.timeOffset);
  const double t2 = static_cast<double>(point2.timeOffset);
  const float v1 = point1.value;
  const float v2 = point2.value;

  const double weight = -t1 / (t2 - t1);
  return v1 + (v2 - v1) * weight;
}

// -----------------------------------------------------------------------------

bool LinearTimeInterpolate::hasGapInRange(
    const std::vector<InterpolationPoint>& candidatePoints,
    const std::vector<InterpolationPoint>& uniquePoints,
    size_t initialUniquePointIdx, size_t subsequentUniquePointIdx) const {
  // Get the time offsets of the interpolation/extrapolation points
  const int64_t t1 = uniquePoints[initialUniquePointIdx].timeOffset;
  const int64_t t2 = uniquePoints[subsequentUniquePointIdx].timeOffset;
  const int64_t targetTime = 0;  // by definition

  // Determine the temporal range including target and both points
  const int64_t minTime = std::min({t1, t2, targetTime});
  const int64_t maxTime = std::max({t1, t2, targetTime});

  // Check if any point in candidatePoints falls within the temporal range and
  // has a missing value
  for (const auto& point : candidatePoints) {
    if (point.timeOffset >= minTime && point.timeOffset <= maxTime &&
        !point.isValid()) {
      return true;  // Gap detected
    }
  }

  return false;  // No gaps found
}

// -----------------------------------------------------------------------------

const ufo::Variables& LinearTimeInterpolate::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
