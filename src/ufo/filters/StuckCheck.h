
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_STUCKCHECK_H_
#define UFO_FILTERS_STUCKCHECK_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/StuckCheckParameters.h"
#include "ufo/filters/TrackCheckUtils.h"

namespace util {
class DateTime;
}

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

/// Flags sequential observations whose filter variables have streaks of unchanging
/// measurements over a time span such that the following two conditions are satisfied:
/// 1. The streak exceeds a user-specified number of observations (\p numberStuckTolerance)
/// 2. The streak continues for longer than a user-specified duration (\p timeStuckTolerance)
/// and/or the streak continues throughout every one of the observations within the observation
/// grouping.
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
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~StuckCheck() override;

 private:
  Parameters_ options_;
  // Instantiate object for accessing observations that may be held on multiple MPI ranks.
  std::unique_ptr<std::vector<util::DateTime>> obsGroupDateTimes_;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::track;}
  std::vector<float> collectStationVariableData(
      std::vector<size_t>::const_iterator stationObsIndicesBegin,
      std::vector<size_t>::const_iterator stationObsIndicesEnd,
      const std::vector<size_t> &validObsIds,
      const std::vector<float> &globalData) const;
  void potentiallyRejectStreak(std::vector<size_t>::const_iterator stationIndicesBegin,
                               std::vector<size_t>::const_iterator stationIndicesEnd,
                               const std::vector<size_t> &validObsIds,
                               size_t startOfStreakIndex,
                               size_t endOfStreakIndex,
                               std::vector<bool> &isRejected,
                               std::string stationId) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_STUCKCHECK_H_
