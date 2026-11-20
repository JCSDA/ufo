/*
 * (C) Copyright 2019 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Range.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsTraits.h"
#include "ufo/operators/timeoper/ObsTimeOperParameters.h"
#include "ufo/operators/timeoper/ObsTimeOperUtil.h"
#include "ufo/SampledLocations.h"

namespace ufo {


//--------------------------------------------------------------------------------------------------
std::vector<std::vector<float>> timeWeightCreate(const ioda::ObsSpace & odb_,
                                    const ObsTimeOperParameters & parameters) {
  util::DateTime windowBegin(odb_.windowStart());
  const util::Duration windowSub = parameters.windowSub.value();
  int64_t windowSubSec = windowSub.toSeconds();

  std::size_t nlocs = odb_.nlocs();

  oops::Log::debug() << "nlocs =    " << nlocs << std::endl;

  std::vector<float> TimeWeightObsAfterState(nlocs, 0.0);

  std::vector<util::DateTime> dateTimeIn(nlocs);
  odb_.get_db("MetaData", "dateTime", dateTimeIn);

  oops::Log::debug() << "dateTime =  " << dateTimeIn[0].toString() << std::endl;

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

  oops::Log::trace() << "TimeOperUtil::timeWeightCreate done" << std::endl;
  return timeWeights;
}
// -----------------------------------------------------------------------------
oops::Locations<ObsTraits> timeOperLocations(const ioda::ObsSpace & odb,
                                             const util::Duration & windowSub) {
  const size_t nlocs = odb.nlocs();
  std::vector<float> lons(nlocs), lats(nlocs);
  std::vector<util::DateTime> times(nlocs);
  odb.get_db("MetaData", "longitude", lons);
  odb.get_db("MetaData", "latitude", lats);
  odb.get_db("MetaData", "dateTime", times);

  const util::DateTime windowBegin = odb.windowStart();
  const int64_t windowSubSec = windowSub.toSeconds();
  std::vector<float> pathLons, pathLats;
  std::vector<util::DateTime> pathTimes;
  std::vector<util::Range<size_t>> pathsGroupedByLocation(nlocs);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const int64_t fromStart = (times[jloc] - windowBegin).toSeconds();
    const util::DateTime before =
      windowBegin + util::Duration((fromStart / windowSubSec) * windowSubSec);
    const util::DateTime after = (before == times[jloc]) ? before : before + windowSub;
    for (const util::DateTime & t : {before, after}) {
      pathLons.push_back(lons[jloc]);
      pathLats.push_back(lats[jloc]);
      pathTimes.push_back(t);
    }
    pathsGroupedByLocation[jloc] = {2 * jloc, 2 * jloc + 2};
  }
  // GetValues ignores paths on an excluded window bound, so move them one second inside
  const std::vector<bool> inWindow = odb.timeWindow().createTimeMask(pathTimes);
  for (size_t jpath = 0; jpath < pathTimes.size(); ++jpath) {
    if (!inWindow[jpath])
      pathTimes[jpath] += util::Duration(pathTimes[jpath] <= windowBegin ? 1 : -1);
  }
  return oops::SampledLocations<ObsTraits>(
        std::make_unique<SampledLocations>(pathLons, pathLats, pathTimes, odb.distribution(),
                                           std::move(pathsGroupedByLocation)));
}

// -----------------------------------------------------------------------------
void timeInterpolate(GeoVaLs & gv, const std::vector<std::vector<float>> & timeWeights) {
  const double missing = util::missingValue<double>();
  const oops::Variables vars = gv.getVars();
  oops::Variables newVars = vars;
  newVars -= gv.getReducedVars();
  std::vector<size_t> nlevs;
  for (const auto & var : newVars) nlevs.push_back(gv.nlevs(var, GeoVaLFormat::SAMPLED));
  gv.addReducedVars(newVars, nlevs);

  std::vector<util::Range<size_t>> paths;
  std::vector<double> before, after;
  for (const auto & var : vars) {
    if (gv.areReducedAndSampledFormatsAliased(var)) continue;  // one path per location
    gv.getProfileIndicesGroupedByLocation(var, paths, GeoVaLFormat::SAMPLED);
    const size_t nlevs = gv.nlevs(var, GeoVaLFormat::SAMPLED);
    for (size_t jloc = 0; jloc < paths.size(); ++jloc) {
      ASSERT(paths[jloc].end - paths[jloc].begin == 2);
      before.assign(nlevs, 0.0);
      after.assign(nlevs, 0.0);
      gv.getProfile(before, var, paths[jloc].begin, GeoVaLFormat::SAMPLED);
      gv.getProfile(after, var, paths[jloc].begin + 1, GeoVaLFormat::SAMPLED);
      for (size_t jlev = 0; jlev < before.size(); ++jlev) {
        before[jlev] = (before[jlev] == missing || after[jlev] == missing) ? missing :
          timeWeights[0][jloc] * before[jlev] + timeWeights[1][jloc] * after[jlev];
      }
      gv.putProfile(before, var, jloc, GeoVaLFormat::REDUCED);
    }
  }
}

// -----------------------------------------------------------------------------
void timeInterpolateAD(GeoVaLs & gv, const GeoVaLs & gvad,
                       const std::vector<std::vector<float>> & timeWeights) {
  const double missing = util::missingValue<double>();
  std::vector<util::Range<size_t>> paths;
  std::vector<double> reduced, sampled;
  for (const auto & var : gv.getVars()) {
    gv.getProfileIndicesGroupedByLocation(var, paths, GeoVaLFormat::SAMPLED);
    const size_t nlevs = gv.nlevs(var, GeoVaLFormat::SAMPLED);
    for (size_t jloc = 0; jloc < paths.size(); ++jloc) {
      reduced.assign(nlevs, 0.0);
      gvad.getProfile(reduced, var, jloc, GeoVaLFormat::REDUCED);
      const size_t npaths = paths[jloc].end - paths[jloc].begin;
      for (size_t jp = 0; jp < npaths; ++jp) {
        const double weight = (npaths == 1) ? 1.0 : timeWeights[jp][jloc];
        sampled.assign(nlevs, 0.0);
        gv.getProfile(sampled, var, paths[jloc].begin + jp, GeoVaLFormat::SAMPLED);
        for (size_t jlev = 0; jlev < sampled.size(); ++jlev) {
          if (reduced[jlev] != missing) sampled[jlev] += weight * reduced[jlev];
        }
        gv.putProfile(sampled, var, paths[jloc].begin + jp, GeoVaLFormat::SAMPLED);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
