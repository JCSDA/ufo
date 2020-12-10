/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKSHIPDIAGNOSTICS_H_
#define UFO_FILTERS_TRACKCHECKSHIPDIAGNOSTICS_H_

#include <utility>
#include <vector>

#include "ufo/filters/TrackCheckShip.h"
namespace ufo {
class TrackCheckShipDiagnostics
{
  typedef std::vector<TrackCheckShip::ObservationStatistics> ObsStatsVec;
  typedef TrackCheckShip::TrackStatistics TrkStats;
  typedef std::pair<ObsStatsVec, TrkStats> SingleTrackInitialCalculationResults;
  typedef std::pair<std::vector<size_t>, int> FirstIterativeRemovalInfo;
  std::vector<SingleTrackInitialCalculationResults> multipleTrackInitialCalculationResults_;
  std::vector<bool> earlyBreaks_;
  std::vector<FirstIterativeRemovalInfo> firstIterativeRemovalInfo_;
  std::vector<double> distanceSum_, distancePrevObsOmitted_, distanceCurrentObsOmitted_,
  timeSum_;  // category nine
  std::vector<double> previousSegmentDistanceProportion_,
  previousObservationDistanceAveragedProportion_, previousSegmentTimeProportion_,
  previousAndFastestSegmentTimeProportion_;  // category ten
 public:
  /// \brief Updates the collection of track diagnostics to include
  /// the calculated values from a new track.
  void storeInitialCalculationResults(SingleTrackInitialCalculationResults
                                      singleTrackInitalCalcResults) {
    multipleTrackInitialCalculationResults_.push_back(singleTrackInitalCalcResults);
  }
  /// \brief Returns the full collection of track diagnostics, separated by
  /// track.
  const std::vector<SingleTrackInitialCalculationResults> &getInitialCalculationResults() const {
    return multipleTrackInitialCalculationResults_;
  }

  /// \brief Stores the indicator as to if the track was deemed not worth
  /// checking after the initial calculations were performed
  void storeEarlyBreakResult(bool result) {
    earlyBreaks_.push_back(result);
  }
  /// \brief Returns the collection of indicators as to which tracks were
  /// deemed not worth checking.
  const std::vector<bool> &getEarlyBreaks() const {
    return earlyBreaks_;
  }

  /// \brief Stores the observation(s) removed on the first iteration of the main removal loop.
  void storeFirstIterativeRemovalInfo(
      const FirstIterativeRemovalInfo &firstIterativeRemovalInfo)
  {
    firstIterativeRemovalInfo_.push_back(firstIterativeRemovalInfo);
  }

  const std::vector<FirstIterativeRemovalInfo> &getFirstIterativeRemovalInfo() const
  {
    return firstIterativeRemovalInfo_;
  }

  std::vector<double> getDistanceSum() const {
    return distanceSum_;
  }

  void storeDistanceSum(const double &distanceSum) {
    distanceSum_.push_back(distanceSum);
  }

  std::vector<double> getDistancePrevObsOmitted() const {
    return distancePrevObsOmitted_;
  }

  void storeDistancePrevObsOmitted(double distancePrevObsOmitted) {
    distancePrevObsOmitted_.push_back(distancePrevObsOmitted);
  }

  std::vector<double> getDistanceCurrentObsOmitted() const {
    return distanceCurrentObsOmitted_;
  }
  void storeDistanceCurrentObsOmitted(const double &distanceCurrentObsOmitted) {
    distanceCurrentObsOmitted_.push_back(distanceCurrentObsOmitted);
  }
  std::vector<double> getTimeSum() const {
    return timeSum_;
  }
  void storeTimeSum(const double &timeSum) {
    timeSum_.push_back(timeSum);
  }
  std::vector<double> getPreviousSegmentDistanceProportion() const {
    return previousSegmentDistanceProportion_;
  }
  void storePreviousSegmentDistanceProportion(const double &previousSegmentDistanceProportion) {
    previousSegmentDistanceProportion_.push_back(previousSegmentDistanceProportion);
  }
  std::vector<double> getPreviousObservationDistanceAveragedProportion() const {
    return previousObservationDistanceAveragedProportion_;
  }
  void storePreviousObservationDistanceAveragedProportion(
      const double &previousObservationDistanceAveragedProportion) {
    previousObservationDistanceAveragedProportion_.push_back(
          previousObservationDistanceAveragedProportion);
  }
  std::vector<double> getPreviousSegmentTimeProportion() const {
    return previousSegmentTimeProportion_;
  }
  void storePreviousSegmentTimeProportion(const double &previousSegmentTimeProportion) {
    previousSegmentTimeProportion_.push_back(previousSegmentTimeProportion);
  }
  std::vector<double> getPreviousAndFastestSegmentTimeProportion() const {
    return previousAndFastestSegmentTimeProportion_;
  }
  void storePreviousAndFastestSegmentTimeProportion(
      const double &previousAndFastestSegmentTimeProportion) {
    previousAndFastestSegmentTimeProportion_.push_back(previousAndFastestSegmentTimeProportion);
  }
};
}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPDIAGNOSTICS_H_
