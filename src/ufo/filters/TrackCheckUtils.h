/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKUTILS_H_
#define UFO_FILTERS_TRACKCHECKUTILS_H_

#include <array>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "ufo/filters/Variable.h"

namespace eckit {
class Configuration;
}

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

inline util::Duration abs(const util::Duration &duration) {
  return duration.toSeconds() >= 0 ? duration : -duration;
}

class RecursiveSplitter;

namespace TrackCheckUtils {
typedef std::array<float, 3> Point;
float distance(const Point &a, const Point &b);

enum class CheckResult : char
{
  FAILED = false,
  PASSED = true,
  SKIPPED
};

enum class SweepResult {NO_MORE_SWEEPS_REQUIRED, ANOTHER_SWEEP_REQUIRED};

/// \brief Locations/times of all observations processed by the track checking filter.
class ObsGroupLocationTimes {
 public:
  std::vector<float> latitudes;
  std::vector<float> longitudes;
  std::vector<util::DateTime> datetimes;
};

class ObsLocationTime {
 public:
  ObsLocationTime(float latitude, float longitude,
                          const util::DateTime &time);
  const Point &location() const { return location_; }
  const util::DateTime &time() const { return time_; }
 private:
  Point location_;
  util::DateTime time_;
};

class CheckCounter {
 public:
  CheckCounter();

  float failedChecksFraction() const;
  void registerCheckResult(const CheckResult &result);
  void unregisterCheckResult(const CheckResult &result);
 private:
  int numFailedChecks_;
  int numChecks_;
};

std::vector<size_t> getValidObservationIds(
    const std::vector<bool> &apply,
    const boost::shared_ptr<ioda::ObsDataVector<int>> &flags);

void groupObservationsByStation(const std::vector<size_t> &validObsIds,
                                RecursiveSplitter &splitter,
                                const eckit::Configuration &config,
                                const ioda::ObsSpace &obsdb);

void groupObservationsByRecordNumber(const std::vector<size_t> &validObsIds,
                                     RecursiveSplitter &splitter,
                                     const ioda::ObsSpace &obsdb);

void groupObservationsByVariable(const Variable &variable,
                                 const std::vector<size_t> &validObsIds,
                                 RecursiveSplitter &splitter, const ioda::ObsSpace &obsdb);

template <typename VariableType>
void groupObservationsByTypedVariable(const Variable &variable,
                                      const std::vector<size_t> &validObsIds,
                                      RecursiveSplitter &splitter,
                                      const ioda::ObsSpace &obsdb);

void sortTracksChronologically(const std::vector<size_t> &validObsIds,
                               RecursiveSplitter &splitter,
                               const ioda::ObsSpace &obsdb);

ObsGroupLocationTimes collectObservationsLocations(const ioda::ObsSpace &obsdb);

void flagRejectedObservations(const std::vector<bool> &isRejected,
                              std::vector<std::vector<bool> > &flagged);

}  // namespace TrackCheckUtils

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKUTILS_H_
