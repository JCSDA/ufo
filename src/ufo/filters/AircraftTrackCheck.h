/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_AIRCRAFTTRACKCHECK_H_
#define UFO_FILTERS_AIRCRAFTTRACKCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
class Configuration;
}

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

class AircraftTrackCheckParameters;
class PiecewiseLinearInterpolation;
class RecursiveSplitter;

/// \brief Checks aircraft tracks, rejecting observations inconsistent with the rest of the track.
///
/// Each track (i.e. the observations taken during each flight) is checked separately. The
/// algorithm performs a series of sweeps over the observations from each track. For each
/// observation, multiple estimates of the instantaneous aircraft speed and ascent/descent rate are
/// obtained by comparing the reported aircraft position with the positions reported during a
/// number a nearby (earlier and later) observations that haven't been rejected in previous sweeps.
/// An observation is rejected if a certain fraction of these estimates lie outside the valid
/// range. Sweeps continue until one of them fails to reject any observations.
///
/// See AircraftTrackCheckParameters for the documentation of this filter's parameters.
class AircraftTrackCheck : public FilterBase,
    private util::ObjectCounter<AircraftTrackCheck> {
public:
  static const std::string classname() {return "ufo::AircraftTrackCheck";}

  AircraftTrackCheck(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                     boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                     boost::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~AircraftTrackCheck() override;

private:
  enum class SweepResult {NO_MORE_SWEEPS_REQUIRED, ANOTHER_SWEEP_REQUIRED};
  struct ObsData;
  struct TrackObservation;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::track;}

  std::vector<size_t> getValidObservationIds(const std::vector<bool> &apply) const;

  void groupObservationsByFlightId(const std::vector<size_t> &validObsIds,
                                   RecursiveSplitter &splitter) const;

  void sortTracksChronologically(const std::vector<size_t> &validObsIds,
                                 RecursiveSplitter &splitter) const;

  ObsData collectObsData() const;

  void checkTracks(const std::vector<size_t> &validObsIds,
                   RecursiveSplitter &splitter) const;

  /// Returns an interpolator mapping pressures (in Pa) to maximum accepted speeds (in km/s).
  PiecewiseLinearInterpolation makeMaxSpeedByPressureInterpolation() const;

  void identifyRejectedObservationsInTrack(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const ObsData &obsData,
      const PiecewiseLinearInterpolation &maxSpeedByPressure,
      std::vector<bool> &isRejected) const;

  std::vector<TrackObservation> collectTrackObservations(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const ObsData &obsData) const;

  /// Iterate once over all observations in \p trackObservations, rejecting those inconsistent
  /// with nearby observations.
  ///
  /// \param[inout] trackObservations
  ///   Attributes of all observations in a track. Modified in place.
  /// \param[in]
  ///   Dependence of the expected maximum aircraft speed on air pressure (and thus height).
  /// \param[inout] workspace
  ///   A vector used internally by the function, passed by parameter to avoid repeated memory
  ///   allocations and deallocations.
  SweepResult sweepOverObservations(
      std::vector<TrackObservation> &trackObservations,
      const PiecewiseLinearInterpolation &maxValidSpeedAtPressure,
      std::vector<float> &workspace) const;

  void flagRejectedTrackObservations(
      std::vector<size_t>::const_iterator trackObsIndicesBegin,
      std::vector<size_t>::const_iterator trackObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const std::vector<TrackObservation> &TrackObservation,
      std::vector<bool> &isRejected) const;

  void flagRejectedObservations(const std::vector<bool> &isRejected,
                                std::vector<std::vector<bool> > &flagged) const;

private:
  std::unique_ptr<AircraftTrackCheckParameters> options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_AIRCRAFTTRACKCHECK_H_
