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
  std::vector<SingleTrackInitialCalculationResults> multipleTrackInitialCalculationResults_;
  std::vector<ObsStatsVec> calculatedResultsSimultaneousDeferred_;
  std::vector<bool> earlyBreaks_;
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

  /// \brief Stores the recalculations of values after deferring simultaneous observations.
  ///
  /// Does not store counter values, because those are not updated after the first iteration.
  void storeCalculatedResultsSimultaneousDeferred(ObsStatsVec obsStatsVec) {
    calculatedResultsSimultaneousDeferred_.push_back(obsStatsVec);
  }

  /// \brief Returns the recalculated values calculated after deferring simultaneous
  /// observations
  const std::vector<ObsStatsVec> &getCalculatedResultsSimultaneousDeferred() const {
    return calculatedResultsSimultaneousDeferred_;
  }
};
}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPDIAGNOSTICS_H_
