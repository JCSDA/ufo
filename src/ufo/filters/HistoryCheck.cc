/*
 * (C) 2021 Crown Copyright Met Office. All rights reserved.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/HistoryCheck.h"

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/none.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>
#include "eckit/config/Configuration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/filters/HistoryCheckParameters.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/StuckCheck.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

HistoryCheck::HistoryCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                            std::shared_ptr<ioda::ObsDataVector<int> > flags,
                            std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::debug() << "HistoryCheck: config = " << options_ << "\n";
}

HistoryCheck::HistoryCheck(ioda::ObsSpace &obsdb,
                           const ufo::HistoryCheck::Parameters_ &parameters,
                           std::shared_ptr<ioda::ObsDataVector<int> > flags,
                           std::shared_ptr<ioda::ObsDataVector<float> > obserr,
                           const eckit::LocalConfiguration &conf)
  : HistoryCheck(obsdb, parameters, flags, obserr)
{
  unitTestConfig_ = conf;
}

/// This filter runs the ship track check and stuck check filters consecutively over an auxiliary
/// obs space (which is assumed to be a superset of \p obsdb_ with an earlier starting time, and
/// possibly a later ending time), before checking which observations have both been (1) flagged
/// by either of the sub-filters from the superset obs space and (2) are located within \p obsdb_
void HistoryCheck::applyFilter(const std::vector<bool> & apply,
                               const Variables & filtervars,
                               std::vector<std::vector<bool> > & flagged) const {
  const util::DateTime widerWindowStart =
    obsdb_.windowStart() - options_.timeBeforeStartOfWindow.value();
  const util::DateTime widerWindowEnd = obsdb_.windowEnd() + options_.timeAfterEndOfWindow.value();
  eckit::LocalConfiguration widerWindowConfig;
  widerWindowConfig.set("begin", widerWindowStart.toString());
  widerWindowConfig.set("end", widerWindowEnd.toString());
  const util::TimeWindow widerTimeWindow(widerWindowConfig);
  // In order to prevent the MPI from distributing the aux spaces's observations to different
  // ranks from the distribution used for obsdb_, widerObsSpace uses the myself communicator
  // for both time and spatial communicators, ensuring that all observations in widerObsSpace
  // are saved to all ranks
  ioda::ObsSpace widerObsSpace(options_.largerObsSpace.value().toConfiguration(),
                               oops::mpi::myself(), widerTimeWindow,
                               oops::mpi::myself());
  if (options_.resetLargerObsSpaceVariables) {  // used for unit testing
    if (unitTestConfig_.has("station_ids_wide")) {
      const std::vector<int> stationIds = unitTestConfig_.getIntVector("station_ids_wide");
      widerObsSpace.put_db("MetaData", "stationIdentification", stationIds);
    } else if (unitTestConfig_.has("station_ids_wide_string")) {
      const std::vector<std::string> stationIds =
          unitTestConfig_.getStringVector("station_ids_wide_string");
      widerObsSpace.put_db("MetaData", "stationIdentification", stationIds);
    }
    if (unitTestConfig_.has("air_temperatures_wide")) {
      const std::vector<float> airTemperatures =
          unitTestConfig_.getFloatVector("air_temperatures_wide");
      widerObsSpace.put_db("ObsValue", "airTemperature", airTemperatures);
    }
  }  // end of manual data entry section used for unit testing
  std::shared_ptr<ioda::ObsDataVector<float>> obserrWide(
        new ioda::ObsDataVector<float>(widerObsSpace, widerObsSpace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflagsWide(
        new ioda::ObsDataVector<int>(widerObsSpace, widerObsSpace.obsvariables()));
  const oops::RequiredParameter<SurfaceObservationSubtype> &subtype =
      options_.surfaceObservationSubtype;
  const boost::optional<TrackCheckShipCoreParameters> &trackOptions =
      options_.trackCheckShipParameters.value();
  // If the observation type is one which the track check ship filter should be run on and
  // the necessary filter parameters were set within the configuration file
  if (subtype != SurfaceObservationSubtype::LNDSYB &&
      subtype != SurfaceObservationSubtype::LNDSYN &&
      trackOptions) {
    // Collecting parameters relevant for running the track check ship filter on the wider obs space
    eckit::LocalConfiguration configTrackCheckShip =
        trackOptions->toConfiguration();
    subtype.serialize(configTrackCheckShip);
    // Importing the base parameters
    options_.TrackCheckUtilsParameters::serialize(configTrackCheckShip);
    ufo::TrackCheckShipParameters trackParams;
    trackParams.deserialize(configTrackCheckShip);
    // Setting up and running the track check ship filter on the wider obs space
    ufo::TrackCheckShip trackCheck(widerObsSpace, trackParams,
                                   qcflagsWide, obserrWide);
    trackCheck.preProcess();
  }
  const boost::optional<StuckCheckCoreParameters> &stuckOptions =
      options_.stuckCheckParameters;
  // If the stuck check filter parameters were set and if the stuck check filter has the potential
  // to flag observations (number of observations is greater than the numberStuckTolerance value)
  if (stuckOptions &&
      widerObsSpace.index().size() >
      stuckOptions->numberStuckTolerance.value().value()) {
    // If the observation subtype is one which the stuck check filter should be run on
    if (subtype != SurfaceObservationSubtype::TEMP &&
        subtype != SurfaceObservationSubtype::BATHY &&
        subtype != SurfaceObservationSubtype::TESAC &&
        subtype != SurfaceObservationSubtype::BUOYPROF) {
      // Collecting the relevant parameters for running the stuck check filter
      eckit::LocalConfiguration configStuckCheck =
          stuckOptions->toConfiguration();
      options_.TrackCheckUtilsParameters::serialize(configStuckCheck);
      ufo::StuckCheckParameters stuckParams;
      stuckParams.deserialize(configStuckCheck);
      // Setting up and running the stuck check filter on the wider obs space
      ufo::StuckCheck stuckCheck(widerObsSpace, stuckParams,
                                 qcflagsWide, obserrWide);
      stuckCheck.preProcess();
    }
  }
  // Creating obs accessors for both obs spaces, assuming the same variable is used for grouping
  // into stations on both obs spaces. 3rd arg: recordsAreSingleObs=false always for History Check.
  ObsAccessor historicalObsAccessor = TrackCheckUtils::createObsAccessor(options_.stationIdVariable,
                                                                   widerObsSpace,
                                                                   false);
  ObsAccessor windowObsAccessor =
      TrackCheckUtils::createObsAccessor(options_.stationIdVariable, obsdb_, false);

  std::vector<util::DateTime> wideDts = historicalObsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "dateTime");
  std::vector<float> wideLats = historicalObsAccessor.getFloatVariableFromObsSpace(
        "MetaData", "latitude");
  std::vector<float> wideLons = historicalObsAccessor.getFloatVariableFromObsSpace(
        "MetaData", "longitude");

  std::map<std::string, int> stationIdMap;

  // Create a "central source of truth" map of integer index values for each
  // string-valued station id across both obs spaces.
  const boost::optional<Variable> &statIdVar = options_.stationIdVariable.value();
  // If the station id variable is in use, and the underlying ids are string-formatted.
  if (statIdVar != boost::none &&
      obsdb_.dtype(statIdVar->group(), statIdVar->variable()) == ioda::ObsDtype::String) {
    // retrieve the station ids from both obs spaces
    std::vector<std::string> wideIds = historicalObsAccessor.getStringVariableFromObsSpace(
          statIdVar->group(), statIdVar->variable());
    std::vector<std::string> windowIds = windowObsAccessor.getStringVariableFromObsSpace(
          statIdVar->group(), statIdVar->variable());
    // append the assimilation window ids to the end of the historical obs space's ids
    std::move(windowIds.begin(), windowIds.end(), std::back_inserter(wideIds));
    // remove all non-unique station ids from the collection of ids
    std::sort(wideIds.begin(), wideIds.end());
    const auto endOfUnique = std::unique(wideIds.begin(), wideIds.end());
    wideIds.erase(endOfUnique, wideIds.end());
    // map all unique station ids across both obs spaces to an integer index
    int index = 0;
    for (const std::string& value : wideIds)
      stationIdMap[value] = index++;
  }

  std::vector<int> wideRecordIds = getStationIds(stationIdMap,
                                                 options_.stationIdVariable.value(),
                                                 widerObsSpace, historicalObsAccessor);

  // obsIdentifierData: all of the MetaData needed to uniquely identify each observation
  // MetaData in use: time stamp, lat/lon coordinates, the station id
  // (actual or integer equivalent), and an additional number for differentiating identical
  // observations
  typedef std::tuple<util::DateTime, float, float, int, size_t> obsIdentifierData;

  // Collect identifiers of each flagged observation from the wider obs space (which was used to
  // run the stuck check and/or the track check ship filters. Increment "differentiator
  // counter" until the full identifier to-be-added is unique within the set.
  std::set<obsIdentifierData> wideFlaggedLocationIds;
  // qc flags are the same across all variables for these filters
  const std::vector<int> &wideFlags = (*qcflagsWide)[0];
  for (size_t i = 0; i < wideFlags.size(); i++) {
    if (wideFlags[i] == QCflags::track) {
      obsIdentifierData obsLabel = {
        wideDts.at(i), wideLats.at(i), wideLons.at(i), wideRecordIds.at(i), 0
      };
      while (wideFlaggedLocationIds.find(obsLabel) != wideFlaggedLocationIds.end()) {
        (std::get<4>(obsLabel))++;
      }
      wideFlaggedLocationIds.insert(obsLabel);
    }
  }
  // Retrieve relevant identifier information for all observations in assimilation window.
  std::vector<util::DateTime> windowDts = windowObsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "dateTime");
  std::vector<float> windowLats = windowObsAccessor.getFloatVariableFromObsSpace(
        "MetaData", "latitude");
  std::vector<float> windowLons = windowObsAccessor.getFloatVariableFromObsSpace(
        "MetaData", "longitude");
  std::vector<int> windowStationIds = getStationIds(stationIdMap,
                                                    options_.stationIdVariable.value(),
                                                    obsdb_, windowObsAccessor);

  // Determine vector of locations at which this filter should be applied.
  // If each independent group of observations is stored entirely on a single MPI rank
  // then this vector will be determined separately for each rank.
  // Otherwise, this vector will be concatenated across all ranks.
  const std::vector <bool> globalApply = windowObsAccessor.getGlobalApply(apply);

  // Map obsIdentifierData for every assimilation window observation to its index within the
  // window's full observation accessor.
  boost::unordered_map<obsIdentifierData, const size_t> locationIdToIndex;
  for (size_t i = 0; i < windowDts.size(); i++) {
    if (!globalApply[i]) continue;
    // Set up all observation labels with differentiator counter initially set to 0
    obsIdentifierData obsLabel = {windowDts.at(i), windowLats.at(i), windowLons.at(i),
                                  windowStationIds.at(i), 0};
    // If the observation label has already been added to the map, increment the differentiator
    // counter until the label is not present within the location id map.
    while (locationIdToIndex.find(obsLabel) != locationIdToIndex.end()) {
      (std::get<4>(obsLabel))++;
    }
    locationIdToIndex.insert(std::pair<obsIdentifierData, const size_t>(obsLabel, i));
  }
  // Comparator defined in order to use set intersection operation for map of indexed
  // assimilated observation ids and set of flagged observations from wider observation spaces
  struct setIntersectionComparator {
    bool operator()(const obsIdentifierData &lhs,
                    const std::pair<const obsIdentifierData, const size_t> &rhs) {
      return lhs < rhs.first;
    }
    bool operator()(const std::pair<const obsIdentifierData, const size_t> &lhs,
                    const obsIdentifierData &rhs) {
      return lhs.first < rhs;
    }
  };

  // Iterate through flagged observations in the historical obs space,
  // finding the observations which are also in the assimilation obs space, and
  // marking the associated indices to flag using the ObsAccessor flagRejectedObservations method.
  // The globalApply vector is used to determine which locations should be flagged
  // based on the where clause.
  std::vector<bool> globalObsToFlag(windowObsAccessor.totalNumObservations(), false);
  for (const obsIdentifierData &id : wideFlaggedLocationIds) {
    if (locationIdToIndex.find(id) != locationIdToIndex.end()) {
      const size_t locToFlag = locationIdToIndex.at(id);
      if (globalApply[locToFlag])
        globalObsToFlag.at(locToFlag) = true;
    }
  }
  windowObsAccessor.flagRejectedObservations(globalObsToFlag, flagged);
}

std::vector<int> HistoryCheck::getStationIds(const std::map<std::string, int> &stringMap,
                                             const boost::optional<Variable> &stationIdVar,
                                             const ioda::ObsSpace &obsdb,
                                             const ObsAccessor &obsacc) const {
  if (stationIdVar == boost::none) {
      if (obsdb.obs_group_vars().empty()) {
        // Observations were not grouped into records.
        // Assume all observations were taken by the same station.
        return std::vector<int>(obsacc.totalNumObservations(), 0);
      } else {
        const std::vector<size_t> &recordNumbers = obsacc.getRecordIds();
        return std::vector<int>(recordNumbers.begin(), recordNumbers.end());
      }
  } else {
    switch (obsdb.dtype(stationIdVar->group(), stationIdVar->variable())) {
    case ioda::ObsDtype::Integer:
    {
      return obsacc.getIntVariableFromObsSpace(stationIdVar->group(),
                                               stationIdVar->variable());
    }

    case ioda::ObsDtype::String:
      {
        std::vector<std::string> stringIds =
            obsacc.getStringVariableFromObsSpace(stationIdVar->group(),
                                                 stationIdVar->variable());
        std::vector<int> ints;
        ints.reserve(stringIds.size());
        std::transform(stringIds.begin(), stringIds.end(), std::back_inserter(ints),
                       [&stringMap](const std::string& value) { return stringMap.at(value); });
        return ints;
      }

    default:
      throw eckit::UserError("Only integer and string variables may be used as station IDs",
                             Here());
    }
  }
}

void HistoryCheck::print(std::ostream & os) const {
  os << "HistoryCheck: config = " << options_ << '\n';
}
}  // namespace ufo
