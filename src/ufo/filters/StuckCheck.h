
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_STUCKCHECK_H_
#define UFO_FILTERS_STUCKCHECK_H_

#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/StuckCheckParameters.h"
#include "ufo/filters/TrackCheckUtils.h"

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

/// Flags sequential observations whose filter variables have streaks of unchanging
/// measurements.
///
/// The filter supports two modes:
/// 1. Count + time mode: a streak is flagged if it exceeds a user-specified number of
///    observations (\p numberStuckTolerance) and continues for longer than a
///    user-specified duration (\p timeStuckTolerance), unless the streak spans the
///    whole record/group.
/// 2. Percentage mode: the number stuck tolerance is derived from
///    \p percentageStuckTolerance for each record.
///
/// Types of observations that this check might apply to include the following:
/// LNDSYN, LNDSYB, SHPSYN, SHPSYB, BUOY, MOBSYN, and OPENROAD
///
class StuckCheck: public FilterBase,
    private util::ObjectCounter<StuckCheck> {
 public:
  typedef StuckCheckParameters Parameters_;

  static const std::string classname() {return "ufo::StuckCheck";}

  StuckCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
             ioda::ObsDataVector<int> & flags,
             ioda::ObsDataVector<float> & obserr);

  ~StuckCheck() override;

 private:
  Parameters_ options_;
  mutable std::vector<util::DateTime> obsGroupDateTimes_;
  mutable std::vector<int> numberStuckToleranceVarValues_;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::track;}
  std::vector<float> collectStationVariableData(
      std::vector<size_t>::const_iterator stationObsIndicesBegin,
      std::vector<size_t>::const_iterator stationObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const std::vector<float> &globalData) const;
  boost::optional<size_t> getStuckToleranceForRecord(
      const std::vector<size_t> &validObsIds,
      std::vector<size_t>::const_iterator stationObsIndicesBegin,
      std::vector<size_t>::const_iterator stationObsIndicesEnd) const;
  void potentiallyRejectStreak(std::vector<size_t>::const_iterator stationIndicesBegin,
                               std::vector<size_t>::const_iterator stationIndicesEnd,
                               const std::vector<size_t> &validObsIds,
                               size_t startOfStreakIndex,
                               size_t endOfStreakIndex,
                               std::vector<bool> &isRejected,
                               std::string stationId,
                               size_t stuckTolerance) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_STUCKCHECK_H_
