/*
 * (C) Copyright 2019 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/timeoper/ObsTimeOperUtil.h"

namespace ufo {


//--------------------------------------------------------------------------------------------------
std::vector<std::vector<float>> timeWeightCreate(const ioda::ObsSpace & odb_,
                                    const eckit::Configuration & config) {
  util::DateTime windowBegin(odb_.windowStart());
  util::Duration windowSub;
  windowSub = util::Duration(config.getString("windowSub"));
  int64_t windowSubSec = windowSub.toSeconds();

  std::size_t nlocs = odb_.nlocs();

  oops::Log::debug() << "nlocs =    " << nlocs << std::endl;

  std::vector<float> TimeWeightObsAfterState(nlocs, 0.0);

  std::vector<util::DateTime> dateTimeIn(nlocs);
  odb_.get_db("MetaData", "datetime", dateTimeIn);

  oops::Log::debug() << "datetime =  " << dateTimeIn[0].toString() << std::endl;

  for (std::size_t i = 0; i < nlocs; ++i) {
    util::Duration timeFromStart = dateTimeIn[i] - windowBegin;
    int64_t timeFromStartSec = timeFromStart.toSeconds();
    int64_t StateTimeFromStartSec =
      (timeFromStartSec / windowSubSec) * windowSubSec;
    if ((timeFromStartSec - StateTimeFromStartSec) == 0) {
      TimeWeightObsAfterState[i] = 1.0f;
    } else {
      TimeWeightObsAfterState[i] = 1.0f - static_cast<float>(timeFromStartSec -
                                                       StateTimeFromStartSec)/
                                         static_cast<float>(windowSubSec);
    }
    oops::Log::debug() << " timeFromStartSec = " << timeFromStartSec
                       << " windowSubSec = " << windowSubSec
                       << " StateTimeFromStartSec = " << StateTimeFromStartSec
                       << std::endl;
  }
  for (std::size_t i=0; i < TimeWeightObsAfterState.size(); ++i) {
    oops::Log::debug() << "timeweights [" << i << "] = "
                       << TimeWeightObsAfterState[i] << std::endl;
  }

  std::vector<float> TimeWeightObsBeforeState(nlocs, 0.0);
  transform(TimeWeightObsAfterState.cbegin(), TimeWeightObsAfterState.cend(),
            TimeWeightObsBeforeState.begin(),
            [] (float element) {return 1.0f - element;});

  std::vector<std::vector<float>> timeWeights;
  timeWeights.push_back(TimeWeightObsAfterState);
  timeWeights.push_back(TimeWeightObsBeforeState);

  for (auto i : timeWeights[0]) {
    oops::Log::debug() << "TimeOperUtil::timeWeights[0] = " << i << std::endl;
  }
  for (auto i : timeWeights[1]) {
    oops::Log::debug() << "TimeOperUtil::timeWeights[1] = " << i << std::endl;
  }

  return timeWeights;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
