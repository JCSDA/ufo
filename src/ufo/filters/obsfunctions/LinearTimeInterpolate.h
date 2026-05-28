/*
 * (C) British Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_LINEARTIMEINTERPOLATE_H_
#define UFO_FILTERS_OBSFUNCTIONS_LINEARTIMEINTERPOLATE_H_

#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling the LinearTimeInterpolate ObsFunction
class LinearTimeInterpolateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearTimeInterpolateParameters, Parameters)

 public:
  /// Input variables for interpolation (must have at least 2 elements)
  oops::RequiredParameter<std::vector<Variable>> inputVariables{
      "input variables",
      "Variables to interpolate between (at least 2 required). Each element "
      "corresponds to the same-indexed element in input timestamps.",
      this};

  /// Input timestamps corresponding to the input variables (must have at least
  /// 2 elements)
  oops::RequiredParameter<std::vector<std::string>> inputTimestamps{
      "input timestamps",
      "Timestamps corresponding to the input variables (same number as input "
      "variables). All entries must be of the same kind: either all datetime "
      "variable references (e.g. MetaData/dateTime) or all literal timestamp "
      "strings in one of the supported formats: YYYY-MM-DDTHH:MM:SSZ (full "
      "datetime) or THH:MM:SSZ (time only). Mixing is not allowed.",
      this};

  /// Target datetime to interpolate at (default: MetaData/dateTime)
  oops::Parameter<Variable> targetDateTime{
      "target datetime",
      "Datetime variable to interpolate to. Defaults to MetaData/dateTime if "
      "not specified.",
      Variable{"MetaData/dateTime"}, this};

  /// Allow interpolation/extrapolation across gaps in data
  oops::Parameter<bool> allowGapInterpolation{
      "allow gap interpolation",
      "If true (default), missing input values are skipped and the nearest "
      "available valid points are used regardless of gaps. When false, "
      "interpolation and extrapolation are disallowed when there are missing "
      "values in the time range being used. When the target time falls "
      "between two timestamps with missing data in between, the output is "
      "set to missing.",
      true, this};
};

// -----------------------------------------------------------------------------

/// \brief Performs piecewise linear time interpolation or extrapolation.
///
/// This function performs piecewise linear interpolation (or extrapolation)
/// based on input values at different times. The interpolation is performed to
/// a target datetime (by default MetaData/dateTime, but can be customized).
///
/// Input variables and timestamps are paired by index: the ith element of input
/// variables corresponds to the ith element of input timestamps.
///
/// Timestamps are supplied via 'input timestamps'. All entries must be of the
/// same kind:
///   - All variable references: e.g. MetaData/dateTime
///   - All literal time-only strings: THH:MM:SSZ (e.g., 'T01:00:00Z')
///   - All literal full datetime strings: YYYY-MM-DDTHH:MM:SSZ
///     (e.g., '2025-12-16T01:10:20Z')
/// Mixing variable references with literal strings is not allowed.
///
/// For time-only literal format, interpolation is performed using only the
/// time-of-day component (hour, minute, second), ignoring date information.
/// This is useful for diurnal patterns where only the time within a 24-hour
/// cycle matters. For full datetime literal format, both date and time are
/// used.
///
/// For two input points, the formula is:
///   result = value1 + (value2 - value1) * (t - t1) / (t2 - t1)
///
/// For multiple input points, piecewise linear interpolation is performed:
///   - Pairs of (value, timestamp) are sorted by timestamp
///   - Find the two timestamps that bracket the target time
///   - Interpolate linearly between those two points
///   - For times before the first timestamp or after the last timestamp,
///     linear extrapolation is used based on the nearest two points
///
/// Missing data handling depends on the 'allow gap interpolation' parameter:
///   - When true (default): Missing input values are skipped.
///   Interpolation/extrapolation
///     uses the nearest available valid points, even if there are gaps in the
///     time series.
///   - When false: Interpolation/extrapolation across gaps is disallowed. If
///     the time range between the two points being used contains any missing
///     values, the output is set to missing.
///
/// Example with missing data: Input variables at times T1, T2, T3, T4 where T3
/// has missing value.
///   - Default mode (true): Target between T2-T3 or T3-T4 interpolates using
///     T2 and T4.
///   - Strict mode (false): Target between T2-T3 or T3-T4 returns missing (gap
///     detected).
///
/// If multiple input values have identical timestamps but different values, the
/// output is set to missing and a warning is logged. Duplicate timestamps with
/// identical values are allowed and treated as a single point.
///
/// Example (2 points):
///   obs function:
///     name: ObsFunction/LinearTimeInterpolate
///     options:
///       input variables:
///         - name: ObsValue/airTemperatureNewer
///         - name: ObsValue/airTemperatureOlder
///       input timestamps:
///         - MetaData/dateTimeNewer
///         - MetaData/dateTimeOlder
///       target datetime: MetaData/dateTime  # optional
///
/// Example (multiple points for piecewise interpolation):
///   obs function:
///     name: ObsFunction/LinearTimeInterpolate
///     options:
///       input variables:
///         - name: ObsValue/temperature00
///         - name: ObsValue/temperature06
///         - name: ObsValue/temperature12
///         - name: ObsValue/temperature18
///       input timestamps:
///         - MetaData/dateTime00
///         - MetaData/dateTime06
///         - MetaData/dateTime12
///         - MetaData/dateTime18
///
/// Example (literal time-only timestamps):
///   obs function:
///     name: ObsFunction/LinearTimeInterpolate
///     options:
///       input variables:
///         - name: ObsValue/airTemperatureT010000
///         - name: ObsValue/airTemperatureT020000
///         - name: ObsValue/airTemperatureT030000
///       input timestamps:
///         - T01:00:00Z
///         - T02:00:00Z
///         - T03:00:00Z
///
/// Example (literal full datetime timestamps):
///   obs function:
///     name: ObsFunction/LinearTimeInterpolate
///     options:
///       input variables:
///         - name: ObsValue/temperature00
///         - name: ObsValue/temperature06
///         - name: ObsValue/temperature12
///       input timestamps:
///         - 2025-12-16T00:00:00Z
///         - 2025-12-16T06:00:00Z
///         - 2025-12-16T12:00:00Z
///
class LinearTimeInterpolate : public ObsFunctionBase<float> {
 public:
  explicit LinearTimeInterpolate(const eckit::LocalConfiguration &);
  ~LinearTimeInterpolate();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables &requiredVariables() const;

 private:
  /// \brief Data structure for a point used in time interpolation
  struct InterpolationPoint {
    int64_t timeOffset;    ///< Time offset from target DateTime in seconds
    float value;           ///< Observation value
    size_t inputVarIndex;  ///< Index of the input variable (0, 1, 2, ...)

    InterpolationPoint(int64_t t, float v, size_t idx)
        : timeOffset(t), value(v), inputVarIndex(idx) {}

    bool isValid() const { return value != util::missingValue<float>(); }
  };

  /// Attempt to parse a supported literal timestamp string into DateTime.
  /// Returns true when parsing succeeds (and outputs
  /// parsedDateTime/isTimeOnly).
  bool tryParseLiteralTimestamp(const std::string &varToken,
                                util::DateTime &parsedDateTime,
                                bool &parsedIsTimeOnly) const;

  /// Check if there are any missing input variables in the range being used for
  /// interpolation. candidatePoints contains all input variables (including
  /// missing),
  /// uniquePoints contains only valid ones.
  bool hasGapInRange(const std::vector<InterpolationPoint> &candidatePoints,
                     const std::vector<InterpolationPoint> &uniquePoints,
                     size_t initialUniquePointIdx,
                     size_t subsequentUniquePointIdx) const;

  /// Perform piecewise linear interpolation using configured timestamp mode
  /// (all literals or all variable references).
  void performInterpolation(
      ioda::ObsDataVector<float> &out,
      const std::vector<util::DateTime> &targetDateTimes, size_t nlocs,
      size_t numPoints,
      const std::vector<ioda::ObsDataVector<float>> &inputVariablesValues,
      const std::vector<ioda::ObsDataVector<util::DateTime>> &inputDateTimes)
      const;

  /// Calculate time offset between input and target timestamps
  int64_t calculateTimeOffset(const util::DateTime &targetDateTime,
                              const util::DateTime &inputDateTime,
                              bool isTimeOnly) const;

  /// Build list of all interpolation points using configured timestamp sources
  /// (either all variable references or all literal timestamps).
  std::vector<InterpolationPoint> buildInterpolationPoints(
      size_t numPoints, size_t ichan, size_t iloc,
      const util::DateTime &targetDateTime,
      const std::vector<ioda::ObsDataVector<float>> &inputVariablesValues,
      const std::vector<ioda::ObsDataVector<util::DateTime>> &inputDateTimes)
      const;

  /// Extract valid (non-missing) points from all points
  std::vector<InterpolationPoint> extractValidPoints(
      const std::vector<InterpolationPoint> &candidatePoints) const;

  /// Remove duplicate timestamps, returning empty vector if conflicts (same
  /// timestamp, different data) detected
  std::vector<InterpolationPoint> deduplicatePoints(
      const std::vector<InterpolationPoint> &validPoints, size_t iloc,
      size_t ichan) const;

  /// Find the indices of the two points that either bracket the target time
  /// or are the nearest two points for extrapolation
  void selectBracketingOrNearestPoints(
      const std::vector<InterpolationPoint> &uniquePoints, size_t &idx1,
      size_t &idx2) const;

  /// Perform linear interpolation or extrapolation between two points at the
  /// target time (which is implicitly at time offset 0)
  float evaluateLinearFitAtTargetTime(const InterpolationPoint &point1,
                                    const InterpolationPoint &point2) const;

  LinearTimeInterpolateParameters options_;
  ufo::Variables invars_;
  bool timestampsAreLiteral_ = false;
  std::vector<util::DateTime> literalTimestamps_;
  bool timestampsAreTimeOnly_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_LINEARTIMEINTERPOLATE_H_
