/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECK_H_
#define UFO_FILTERS_TRACKCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
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

class ObsAccessor;
class PiecewiseLinearInterpolation;
class RecursiveSplitter;
class TrackCheckParameters;

/// \brief Checks tracks of mobile weather stations, rejecting observations inconsistent with the
/// rest of the track.
///
/// Each track is checked separately. The algorithm performs a series of sweeps over the
/// observations from each track. For each observation, multiple estimates of the instantaneous
/// speed and ascent/descent rate are obtained by comparing the reported position with the
/// positions reported during a number a nearby (earlier and later) observations that haven't been
/// rejected in previous sweeps. An observation is rejected if a certain fraction of these
/// estimates lie outside the valid range. Sweeps continue until one of them fails to reject any
/// observations, i.e. the set of retained observations is self-consistent.
///
/// See TrackCheckParameters for the documentation of this filter's parameters.
///
/// Note: this filter was originally written with aircraft observations in mind. However, it can
/// potentially be useful also for other observation types.
///
class TrackCheck : public FilterBase,
    private util::ObjectCounter<TrackCheck> {
  enum Direction { FORWARD, BACKWARD, NUM_DIRECTIONS };

  /// \brief Results of cross-checking an observation with another (a "buddy").
  struct CheckResults {
    CheckResults() :
      isBuddyDistinct(false),
      speedCheckResult(TrackCheckUtils::CheckResult::SKIPPED),
      climbRateCheckResult(TrackCheckUtils::CheckResult::SKIPPED) {}

    bool isBuddyDistinct;
    TrackCheckUtils::CheckResult speedCheckResult;
    TrackCheckUtils::CheckResult climbRateCheckResult;
  };

  static const int NO_PREVIOUS_SWEEP = -1;

  struct ObsGroupPressureLocationTime {
    TrackCheckUtils::ObsGroupLocationTimes locationTimes;
    std::vector<float> pressures;
  };

  ObsGroupPressureLocationTime collectObsPressuresLocationsTimes(
      const ObsAccessor &obsAccessor) const;

 public:
  static const std::string classname() { return "ufo::TrackCheck"; }

  TrackCheck(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
             std::shared_ptr<ioda::ObsDataVector<int> > flags,
             std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~TrackCheck() override;

 private:
  /// \brief Attributes of an observation belonging to a track.
  class TrackObservation {
   public:
    TrackObservation(float latitude, float longitude, const util::DateTime &time, float pressure);
    float pressure() const { return pressure_; }
    bool rejectedInPreviousSweep() const { return rejectedInPreviousSweep_; }
    bool rejectedBeforePreviousSweep() const { return rejectedBeforePreviousSweep_; }
    bool rejected() const {
      return rejectedInPreviousSweep_ || rejectedBeforePreviousSweep_;
    }
    int numNeighborsVisitedInPreviousSweep(Direction dir) const {
      return numNeighborsVisitedInPreviousSweep_[dir];
    }
    void setNumNeighborsVisitedInPreviousSweep(Direction dir, int num) {
      numNeighborsVisitedInPreviousSweep_[dir] = num;
    }

    /// Estimates the instantaneous speed and climb rate by comparing this observation against
    /// \p buddyObs.
    /// Checks if these estimates are in the accepted ranges and if the two observations
    /// are far enough from each other to be considered "distinct".
    ///
    /// \param buddyObs Observation to compare against.
    /// \param options Track check options.
    /// \param maxValidSpeedAtPressure
    ///   Function mapping air pressure (in Pa) to the maximum realistic speed (in m/s).
    /// \param referencePressure
    ///   Pressure at which the maximum speed should be evaluated.
    ///
    /// \returns An object enapsulating the check results.
    CheckResults checkAgainstBuddy(const TrackObservation &buddyObs,
                                   const TrackCheckParameters &options,
                                   const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
                                   float referencePressure) const;
    void registerCheckResults(const CheckResults &result);
    void unregisterCheckResults(const CheckResults &result);
    void registerSweepOutcome(bool rejectedInSweep);
    float getFailedChecksFraction();

   private:
    TrackCheckUtils::ObsLocationTime obsLocationTime_;
    TrackCheckUtils::CheckCounter checkCounter_;
    float pressure_;
    bool rejectedInPreviousSweep_;
    bool rejectedBeforePreviousSweep_;
    int numNeighborsVisitedInPreviousSweep_[NUM_DIRECTIONS];
  };

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

  ObsAccessor createObsAccessor() const;

  /// Returns an interpolator mapping pressures (in Pa) to maximum accepted speeds (in km/s).
  PiecewiseLinearInterpolation makeMaxSpeedByPressureInterpolation() const;

  void identifyRejectedObservationsInTrack(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const ObsGroupPressureLocationTime &obsPressureLoc,
      const PiecewiseLinearInterpolation &maxSpeedByPressure,
      std::vector<bool> &isRejected) const;

  std::vector<TrackObservation> collectTrackObservations(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const ObsGroupPressureLocationTime &obsPressureLoc) const;

  /// Iterate once over all observations in \p trackObservations, rejecting those inconsistent
  /// with nearby observations.
  ///
  /// \param[inout] trackObservations
  ///   Attributes of all observations in a track. Modified in place.
  /// \param[in]
  ///   Dependence of the expected maximum speed on air pressure (and thus height).
  /// \param[inout] workspace
  ///   A vector used internally by the function, passed by parameter to avoid repeated memory
  ///   allocations and deallocations.
  TrackCheckUtils::SweepResult sweepOverObservations(
      std::vector<TrackObservation> &trackObservations,
      const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
      std::vector<float> &workspace) const;

 private:
  std::unique_ptr<TrackCheckParameters> options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECK_H_
