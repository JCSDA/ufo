/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeDuplicateCheck.h"

#include <algorithm>

#include "ufo/filters/ObsAccessor.h"
#include "ufo/utils/metoffice/MetOfficeSort.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

MetOfficeDuplicateCheck::MetOfficeDuplicateCheck
(ioda::ObsSpace & obsdb,
 const Parameters_ &parameters,
 std::shared_ptr<ioda::ObsDataVector<int> > flags,
 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::debug() << "MetOfficeDuplicateCheck: config = " << options_ << std::endl;
}

// Required for the correct destruction of options_.
MetOfficeDuplicateCheck::~MetOfficeDuplicateCheck()
{}

void MetOfficeDuplicateCheck::applyFilter(const std::vector<bool> & apply,
                                          const Variables & filtervars,
                                          std::vector<std::vector<bool>> & flagged) const {
  ObsAccessor obsAccessor = createObsAccessor();

  const int nlocs = obsAccessor.totalNumObservations();

  std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply);
  if (options_.preSortByRecordID.value()) {
    // The user wants to process observations in fixed (non-random) order. Ensure the filter
    // produces the same results regardless of the number of MPI ranks by ordering the observations
    // to be processed as if we were running in serial: by record ID.
    const std::vector<size_t> recordIds = obsAccessor.getRecordIds();
    std::stable_sort(validObsIds.begin(), validObsIds.end(),
                     [&recordIds](size_t obsIdA, size_t obsIdB)
                     { return recordIds[obsIdA] < recordIds[obsIdB]; });
  }
  std::vector<bool> isThinned(nlocs, false);

  const std::vector <float> latitude =
    obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
  const std::vector <float> longitude =
    obsAccessor.getFloatVariableFromObsSpace("MetaData", "longitude");
  const std::vector <util::DateTime> dateTime =
    obsAccessor.getDateTimeVariableFromObsSpace("MetaData", "dateTime");
  // Convert dateTimes to offsets relative to start of time window.
  std::vector <int64_t> dateTimeOffset(dateTime.size());
  const util::DateTime & winstart = obsdb_.windowStart();
  std::transform(dateTime.cbegin(), dateTime.cend(),
                 dateTimeOffset.begin(),
                 [&winstart](util::DateTime d){return (d - winstart).toSeconds();});
  Variable pressureVariable = options_.verticalCoordinateVariableName.value();
  const std::vector <float> pressure =
    obsAccessor.getFloatVariableFromObsSpace(pressureVariable.group(), pressureVariable.variable());
  std::vector <int> priority =
    obsAccessor.getIntVariableFromObsSpace("MetaData", options_.priorityName);

  // Optional category variable. If defined, determine the locations associated with
  // each distinct value of this variable.
  std::vector <std::string> categoryVariable;
  std::unordered_map <std::string, std::vector<int>> categoryVariableLocations;
  if (options_.categoryVariableName.value() != boost::none) {
    categoryVariable = obsAccessor.getStringVariableFromObsSpace
      ("MetaData", options_.categoryVariableName.value().value());
    for (size_t obsId : validObsIds)
      categoryVariableLocations[categoryVariable[obsId]].push_back(obsId);
  }

  const float latBandWidth = options_.latBandWidth.value();
  const int numLatitudeBands = std::ceil(180.0 / latBandWidth);
  // Assign a latitude band index to each location.
  std::vector <int> latitudeBand;
  for (int jlat = 0; jlat < latitude.size(); ++jlat)
    latitudeBand.push_back(std::max(0, static_cast<int>((89.999 - latitude[jlat]) / latBandWidth)));

  // Sort on latitude band.
  metOfficeSort(validObsIds.begin(), validObsIds.end(),
                [&latitudeBand] (size_t id) {return latitudeBand[id];});
  std::vector <int> latitudeBandSorted;
  for (int jloc = 0; jloc < validObsIds.size(); ++jloc)
    latitudeBandSorted.push_back(latitudeBand[validObsIds[jloc]]);

  // Split into groups of sorted latitude band
  // (the second argument enables OPS compatibility mode).
  RecursiveSplitter splitter(validObsIds.size(), true);
  splitter.groupBy(latitudeBandSorted);
  // Final vector of locations sorted by latitude band and then longitude.
  std::vector<size_t> sortedFinal;
  for (auto categoryGroup : splitter.groups()) {
    std::vector<size_t> obsIdsInCategory;
    for (size_t validObsIndex : categoryGroup)
      obsIdsInCategory.push_back(validObsIds[validObsIndex]);
    RecursiveSplitter catsplitter(obsIdsInCategory.size());
    // Sort by longitude in each latitude band group.
    catsplitter.sortGroupsBy([&longitude, &obsIdsInCategory](size_t ind)
                             {return longitude[obsIdsInCategory[ind]];});
    std::vector<size_t> obsIdsSorted;
    for (auto catgroup : catsplitter.groups()) {
      for (size_t obsIndex : catgroup) {
        const size_t obsId = obsIdsInCategory[obsIndex];
        obsIdsSorted.push_back(obsId);
      }
    }
    sortedFinal.insert(sortedFinal.end(), obsIdsSorted.cbegin(), obsIdsSorted.cend());
  }

  // Location of the start of each latitude band, accounting for empty bands.
  std::vector<int> latBandStartIndex(numLatitudeBands + 1);
  int lband = 0;
  for (int jloc = 0; jloc < sortedFinal.size(); ++jloc) {
    const int jsort = sortedFinal[jloc];
    for (int jband = lband; jband < latitudeBand[jsort] + 1; ++jband)
      latBandStartIndex[jband] = jloc;
    lband = latitudeBand[jsort] + 1;
  }
  for (int jband = lband; jband < numLatitudeBands + 1; ++jband)
    latBandStartIndex[jband] = sortedFinal.size();

  // Spatial and temporal extent of volume used to search for duplicates.
  const float latBinHalfWidth = options_.latBinHalfWidth.value();
  const float lonBinHalfWidth = options_.lonBinHalfWidth.value();
  const float timeBinHalfWidth = options_.timeBinHalfWidth.value();
  const float pBinHalfWidth = options_.pBinHalfWidth.value();

  // Loop over latitude bands and search for pairs of duplicates, rejecting one
  // observation in each pair.
  for (int jband = 0; jband < numLatitudeBands; ++jband) {
    // Vector of first location in each band. Reset for iteration of this loop.
    std::vector<int> latBandCurrentStart(latBandStartIndex);
    // Loop through locations in this band.
    for (int jloc = latBandStartIndex[jband]; jloc < latBandStartIndex[jband + 1]; ++jloc) {
      const int jlocSort = sortedFinal[jloc];
      if (priority[jlocSort] <= 0) continue;  // Move to next location in this band.
      // Search for duplicates in both the current band and the next one (if it is present).
      for (int jbandSearch = jband;
           jbandSearch <= jband + 1 && jbandSearch <= numLatitudeBands; ++jbandSearch) {
        // This ensures the loop over bands can be exited correctly.
        bool nextjloc = false;
        // Indicates that this is the first location in the search band.
        bool firstLocInBand = true;
        // Location at which to start the search.
        const int jlocStart = jbandSearch == jband ? jloc + 1 : latBandCurrentStart[jbandSearch];
        // Loop through locations in the current search band.
        for (int jlocSearch = jlocStart;
             jlocSearch < latBandStartIndex[jbandSearch + 1]; ++jlocSearch) {
          const int jlocSearchSort = sortedFinal[jlocSearch];
          // Lower and upper longitude bounds for this location.
          const float lonLowerBound = longitude[jlocSort] - lonBinHalfWidth;
          const float lonUpperBound = longitude[jlocSort] + lonBinHalfWidth;
          if (longitude[jlocSearchSort] > lonUpperBound) break;  // Move to next search band.
          if (longitude[jlocSearchSort] < lonLowerBound) continue;  // Move to next location.
          if (priority[jlocSearchSort] <= 0) continue;  // Move to next location.
          // Potentially update the starting index for the current search band.
          // This improves the efficiency of the algorithm.
          if (firstLocInBand) {
            latBandCurrentStart[jbandSearch] = jlocSearch;
            firstLocInBand = false;
          }
          // Search for duplicates inside the defined volume.
          if (std::fabs(latitude[jlocSort] - latitude[jlocSearchSort]) <= latBinHalfWidth &&
              std::fabs(longitude[jlocSort] - longitude[jlocSearchSort]) <= lonBinHalfWidth &&
              std::fabs(dateTimeOffset[jlocSort] - dateTimeOffset[jlocSearchSort]) <=
              timeBinHalfWidth &&
              std::fabs(pressure[jlocSort] - pressure[jlocSearchSort]) <= pBinHalfWidth) {
            // A duplicate pair of observations has been found. Retain one observation based on
            // values of priority.
            const size_t jlocThin = priority[jlocSort] < priority[jlocSearchSort] ?
                                                         jlocSort : jlocSearchSort;
            // If the category variable has been defined, thin all observations whose category
            // match that of the initial thinned observation.
            if (options_.categoryVariableName.value() != boost::none) {
              const std::string thinnedCategory = categoryVariable[jlocThin];
              for (size_t thinnedLocation : categoryVariableLocations[thinnedCategory]) {
                priority[thinnedLocation] = -1;
                isThinned[thinnedLocation] = true;
              }
            } else {
              priority[jlocThin] = -1;
              isThinned[jlocThin] = true;
            }
            if (priority[jlocSort] < priority[jlocSearchSort]) {
              nextjloc = true;
              break;
            }
          }
        }
        // Move to the next location in the original band.
        if (nextjloc) break;
      }
    }
  }

  // Flag any observations that have been thinned.
  obsAccessor.flagRejectedObservations(isThinned, flagged);
}


ObsAccessor MetOfficeDuplicateCheck::createObsAccessor() const {
  return ObsAccessor::toAllObservations(obsdb_);
}

void MetOfficeDuplicateCheck::print(std::ostream & os) const {
  os << "MetOfficeDuplicateCheck: config = " << options_ << std::endl;
}

}  // namespace ufo
