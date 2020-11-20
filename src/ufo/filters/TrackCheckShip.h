
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKSHIP_H_
#define UFO_FILTERS_TRACKCHECKSHIP_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <list>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/TrackCheckUtils.h"

namespace eckit {
class Configuration;
}

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

class TrackCheckShipDiagnostics;
class TrackCheckShipParameters;
class RecursiveSplitter;
/// \brief Checks tracks of ships and buoys, rejecting observations whose locations
/// and timestamps make them inconsistent with the rest of the track. The full track is rejected
/// if there are too many inconsistent observations.
///
/// Each track is checked separately. The algorithm will first calculate speeds and distances
/// between every two adjacent and alternating observations, and angles between any three adjacent
/// observations. Based on these values, it will increment a set of "error counters" that reflect
/// how many errors exist within the track. If the "early break check" parameter is set to true and
/// the "error counters" sum up to a value greater than or equal to half of the segments, the
/// filter will ignore the current track and move on to the next one. Otherwise, observations will
/// be iteratively removed based on the location of the fastest segment until all segments between
/// observations are slower than the speed set by the "max speed" parameter (or 80% of it if angles
/// between segments are large).
///
/// See TrackCheckShipsParameters for the documentation of this filter's parameters.
///
class TrackCheckShip: public FilterBase,
    private util::ObjectCounter<TrackCheckShip> {
 public:
  static const std::string classname() {return "ufo::TrackCheckShip";}

  TrackCheckShip(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~TrackCheckShip() override;

  /// \brief Relevant calculated values that are specific to each pair or triplet
  /// of adjacent/alternating observations
  struct ObservationStatistics {
    /// \brief \p timeDifference between the same-index observation and the
    /// previous one.
    util::Duration timeDifference{};
    /// \p distance between the same-index observation and the
    /// previous-index observation.
    double distance{};
    /// \p speed between the same-index observation and the
    /// previous-index observation
    double speed{};
    /// \p angle formed by the track between
    /// the previous-index observation, same-index observation,
    /// and next-index observation.
    double angle{};
    /// \p speedAveraged is the calculated speed
    /// between the previous-index observation and the next-index observation (ignoring
    /// the same-index observation).
    double speedAveraged{};
    /// \p distanceAveraged is the calculated distance between the previous-index observation
    /// and the next-index observation (ignoring the same-index observation).
    double distanceAveraged{};
  };

  /// \brief A container for all track-wise counters and calculations
  /// that indicate the overall quality of the tracks' data.
  struct TrackStatistics {
    /// \p numRejectedObs_ is currently unused, but will be updated in
    /// later pull requests.
    int numRejectedObs_{};
    /// \p numShort_ is the number of observations within an hour of the
    /// previous observation.
    int numShort_{};
    /// \p numFast_ is the number of observations with a speed exceeding
    /// a maximum value, which is to be set in the TrackCheckShipParameters
    /// file.
    int numFast_{};
    /// \p numBends_ is a count of 3-consecutive-observation track segments
    /// in which the track shows a 90-degree or greater turn.
    int numBends_{};
    /// \p sumSpeed_ is the sum of all speed values recorded for observations
    /// that have neither been deemed "fast" or "short"
    double sumSpeed_{};
    /// \p meanSpeed_ is the mean of all speed values included in \p sumSpeed_
    double meanSpeed_{};
  };

 private:
  class TrackObservation {
   public:
    TrackObservation(double latitude, double longitude,
                     const util::DateTime &time,
                     const std::shared_ptr<TrackStatistics> & trackStatistics,
                     const std::shared_ptr<TrackCheckUtils::CheckCounter> & checkCounter,
                     size_t observationNumber);
    const TrackCheckUtils::Point& getLocation() const {
      return obsLocationTime_.location();
    }
    const util::DateTime& getTime() const {
      return obsLocationTime_.time();
    }
    void setDistance(double dist);
    void setTimeDifference(util::Duration tDiff);
    void setSpeed(double speed);
    void setAngle(double angle);
    void setDistanceAveraged(double distAvg);
    void setSpeedAveraged(double speedAvg);

    void calculateTwoObservationValues(
        TrackObservation& prevObs, bool firstIteration,
        const TrackCheckShipParameters& options);
    void resetObservationCalculations();
    void calculateThreeObservationValues(
        const TrackObservation& prevObs, const TrackObservation& nextObs,
        bool firstIteration, const TrackCheckShipParameters& options);
    void adjustTwoObservationStatistics(const TrackCheckShipParameters& options) const;
    void adjustThreeObservationStatistics() const;

    const std::shared_ptr<TrackStatistics> getFullTrackStatistics() const;

    const ObservationStatistics &getObservationStatistics() const;

    void setRejected(bool rejectedInSweep) {
      rejected_ = rejectedInSweep;
    }

    size_t getObservationNumber() const;

    bool rejected() const {
      return rejected_;
    }

   private:
    std::shared_ptr<TrackStatistics> fullTrackStatistics_;
    ObservationStatistics observationStatistics_;
    TrackCheckUtils::ObsLocationTime obsLocationTime_;
    size_t observationNumber_;
    bool rejected_;
  };

 public:
  /// \brief Extends \p distance for \p TrackObservation inputs \p a and \p b.
  static double distance(const TrackObservation &a, const TrackObservation &b) {
    return TrackCheckUtils::distance(a.getLocation(), b.getLocation());
  }

  static double speedEstimate(
      const TrackCheckShip::TrackObservation &obs1,
      const TrackCheckShip::TrackObservation &obs2,
      const TrackCheckShipParameters& options);

  static float angle(const TrackCheckShip::TrackObservation &a,
                     const TrackCheckShip::TrackObservation &b,
                     const TrackCheckShip::TrackObservation &c);
  const TrackCheckShipDiagnostics* diagnostics() const;

 private:
  std::unique_ptr<TrackCheckShipParameters> options_;
  std::unique_ptr<TrackCheckShipDiagnostics> diagnostics_;

  void flagRejectedTrackObservations(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const std::vector<TrackObservation> &trackObservations,
      std::vector<bool> &isRejected) const;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::track;}

  enum CalculationMethod { FIRSTITERATION, MAINLOOP };

  void calculateTrackSegmentProperties(
      const std::vector<std::reference_wrapper<TrackObservation>> &trackObservations,
      CalculationMethod calculationMethod = MAINLOOP) const;

  std::vector<TrackObservation> collectTrackObservations(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const TrackCheckUtils::ObsGroupLocationTimes &obsLoc) const;

  bool earlyBreak(const std::vector<std::reference_wrapper<TrackObservation>> &trackObs,
                  const std::string trackId) const;

  void removeFaultyObservation(
      std::vector<std::reference_wrapper<TrackObservation> > &track,
      const std::vector<std::reference_wrapper<TrackObservation> >::iterator &it,
      bool firstIterativeRemoval, const std::string trackId) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKSHIP_H_
