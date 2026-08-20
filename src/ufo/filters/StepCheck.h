/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_STEPCHECK_H_
#define UFO_FILTERS_STEPCHECK_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/StepCheckParameters.h"
#include "ufo/filters/TrackCheckUtils.h"

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
}  // namespace ioda

namespace ufo {

/// Flags observations when there are large jumps ("steps") between successive
/// data points. A step is detected when the absolute difference between
/// consecutive observations exceeds a specified threshold.
///
/// Features:
/// - Basic step flagging with inclusive/exclusive threshold comparison
/// - Number or percentage-based step tolerance (minimum steps before flagging)
/// - Circular difference calculation for periodic data (e.g., wind direction)
/// - Chunking: split data into chunks, calculate mean step within each chunk,
///   then average across chunks to compare against threshold
/// - Stuck chunk removal: ignore chunks where all values are identical
/// - Average step mode: compare average step magnitude to threshold
///
/// Types of observations that this check might apply to include:
/// Surface observations with temporal continuity such as buoy, aircraft, and
/// ship data.
///
class StepCheck : public FilterBase, private util::ObjectCounter<StepCheck> {
 public:
  typedef StepCheckParameters Parameters_;

  static const std::string classname() { return "ufo::StepCheck"; }

  StepCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
            ioda::ObsDataVector<int> & flags,
            ioda::ObsDataVector<float> & obserr);

  ~StepCheck() override;

 private:
  Parameters_ options_;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::step; }

  /// Collect variable data for a single station/record
  std::vector<float> collectStationVariableData(
      std::vector<size_t>::const_iterator stationObsIndicesBegin,
      std::vector<size_t>::const_iterator stationObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const std::vector<float> &globalData) const;

  /// Calculate the difference between two values, optionally using circular
  /// difference
  float calculateDifference(const float val1, const float val2) const;

  /// Calculate the step (difference) between each consecutive pair of values
  std::vector<float> calculateConsecutiveSteps(
      const std::vector<float> &values) const;

  /// Check if a value exceeds the step threshold (using inclusive/exclusive
  /// comparison as configured)
  bool exceedsThreshold(const float value, const float threshold) const;

  /// Process a station's data to detect steps
  void processStation(std::vector<size_t>::const_iterator stationBegin,
                      std::vector<size_t>::const_iterator stationEnd,
                      const std::vector<size_t> &validObsIds,
                      const std::vector<float> &variableData,
                      std::vector<bool> &isRejected) const;

  /// Calculate the average step for a station's data. Called only when
  /// useAverageStep is true. If chunking is enabled, splits data into chunks,
  /// calculates mean step per chunk (removing stuck chunks if requested), then
  /// averages across chunks. Returns -1.0 if there is insufficient data.
  float calculateAverageStep(const std::vector<float> &validData) const;

  /// Determine which observations should be flagged based on threshold and
  /// tolerance logic. For average step mode, calculates the average step
  /// (with chunking if configured) and returns all indices if it exceeds the
  /// threshold. Otherwise, returns indices of observations where individual
  /// steps exceed threshold (subject to tolerance).
  std::vector<size_t> determineObservationsToFlag(
      const std::vector<float> &values) const;

  /// Determine whether to flag observations based on tolerance settings
  bool shouldFlagBasedOnTolerance(const size_t exceedingCount,
                                  const size_t totalSteps) const;

  /// Check if a chunk should be considered "stuck" (all values the same within
  /// tolerance)
  bool isChunkStuck(const std::vector<float> &chunkValues) const;

  /// Calculate the mean of a vector of values
  float calculateMean(const std::vector<float> &values) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_STEPCHECK_H_
