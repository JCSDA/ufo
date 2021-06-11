/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKUTILS_H_
#define UFO_FILTERS_TRACKCHECKUTILS_H_

#include <array>
#include <memory>
#include <vector>

#include <boost/optional.hpp>

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

class ObsAccessor;
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

/// \brief Create an ObsAccessor object providing access to observations that need to be checked
/// by the current MPI task.
///
/// The composition of this set of observations depends on the value of \p stationIdVariable.
///
/// If \p stationIdVariable is empty, observations are assumed to be grouped into tracks by the
/// record ID. Each MPI rank is guaranteed to hold either all or no observations from a given
/// record. Thus the returned ObsAccessor gives access to all observations from the records held on
/// the current MPI rank.
///
/// Otherwise, observations are assumed to be grouped into tracks by the variable \p
/// *stationIdVariable.  If this variable was also used to group observations into records, the
/// returned ObsAccessor is constructed as if \p stationIdVariable was empty; otherwise, it gives
/// access to observations held on all MPI ranks.
ObsAccessor createObsAccessor(const boost::optional<Variable> &stationIdVariable,
                              const ioda::ObsSpace &obsdb);

void sortTracksChronologically(const std::vector<size_t> &validObsIds,
                               const ObsAccessor &obsAccessor,
                               RecursiveSplitter &splitter);

ObsGroupLocationTimes collectObservationsLocations(const ObsAccessor &obsAccessor);

}  // namespace TrackCheckUtils

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKUTILS_H_
