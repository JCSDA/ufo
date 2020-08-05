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
  typedef std::pair<ObsStatsVec, TrkStats> SingleTrackDiagnostics;
  std::vector<SingleTrackDiagnostics> fullRunDiagnostics_;
 public:
  /// \brief Updates the collection of track diagnostics to include
  /// the calculated values from a new track.
  void storeDiagnostics(SingleTrackDiagnostics singleTrackDiagnostics) {
    fullRunDiagnostics_.push_back(singleTrackDiagnostics);
  }
  /// \brief Returns the full collection of track diagnostics, separated by
  /// track.
  const std::vector<SingleTrackDiagnostics> &getDiagnostics() const {
    return fullRunDiagnostics_;
  }
};
}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIPDIAGNOSTICS_H_
