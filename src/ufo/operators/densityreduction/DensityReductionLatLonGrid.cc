/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/densityreduction/DensityReductionLatLonGrid.h"

#include <unordered_map>

#include "oops/mpi/mpi.h"

namespace ufo {

static DensityReductionMaker<DensityReductionLatLonGrid>
makerDensityReductionLatLonGrid_("lat lon grid");

DensityReductionLatLonGrid::DensityReductionLatLonGrid(const Parameters_ & params,
                                                       const ioda::ObsSpace & odb)
  : DensityReductionBase(params, odb),
    params_(params)
{
  // Parameter validation
  if (params_.latMax.value() <= params_.latMin.value()) {
    throw eckit::UserError("Latitude upper bound must be greater than lower bound", Here());
  }
}

void DensityReductionLatLonGrid::fillModifiedLocations(std::vector<float> & latsModified,
                                                       std::vector<float> & lonsModified,
                                                       std::vector<util::DateTime> & timesModified,
                                                       std::vector<util::Range<size_t>> &
                                                       pathsGroupedByLocation) const {
  const int nlocs = odb_.nlocs();
  const int gnlocs = odb_.globalNumLocs();
  const int rank = static_cast<int>(odb_.comm().rank());

  const float dLat = params_.dLat.value();
  const float dLon = params_.dLon.value();
  const float latMin = params_.latMin.value();
  const float latMax = params_.latMax.value();
  const float lonMin = params_.lonMin.value();

  // Compute number of latitude bins.
  const float latSpan = latMax - latMin;
  // Adding one accounts for observations located exactly on the upper boundary.
  const int numBinsLat = static_cast<int>(latSpan / dLat) + 1;

  // Get lat/lon/time from ObsSpace.
  std::vector<float> lats(nlocs);
  std::vector<float> lons(nlocs);
  std::vector<util::DateTime> times(nlocs);
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "dateTime", times);

  // Vector signifying whether each location on this rank is a 'patch' observation.
  // The contents of this vector depends on the ObsSpace distribution used.
  std::vector<bool> patchObsVec(nlocs);
  odb_.distribution()->patchObs(patchObsVec);

  // Fill vector of ranks of patch observations.
  // The contents of this vector depends on the ObsSpace distribution used.
  std::vector<int> patchRanks;
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (patchObsVec[jloc]) {
      patchRanks.push_back(rank);
    }
  }

  // Gather ranks of patch observations.
  // The use of oops::mpi::allGatherv() ensures the gather procedure
  // works for subclasses of ioda::Distribution whose implementation of allGatherv() may
  // not produce the desired result.
  std::vector<int> globalPatchRanks(patchRanks);
  oops::mpi::allGatherv(odb_.comm(), globalPatchRanks);
  ASSERT(globalPatchRanks.size() == gnlocs);

  // Produce global vectors of lat/lon/time
  std::vector<float> globalLatitude(lats);
  odb_.distribution()->allGatherv(globalLatitude);
  std::vector<float> globalLongitude(lons);
  odb_.distribution()->allGatherv(globalLongitude);
  std::vector<util::DateTime> globalTime(times);
  odb_.distribution()->allGatherv(globalTime);

  // Assign each location a bin index in the lat-lon grid.
  // The grid is flattened to 1D.
  std::vector<int64_t> bin(gnlocs);
  for (int gloc = 0; gloc < gnlocs; ++gloc) {
    const int64_t binLat =
      static_cast<int64_t>(std::trunc((globalLatitude[gloc] - latMin) / dLat));
    const int64_t binLon =
      static_cast<int64_t>(std::trunc((globalLongitude[gloc] - lonMin) / dLon));
    bin[gloc] = binLat + numBinsLat * binLon;
  }

  // Determine all locations associated with each populated lat-lon bin.
  // Unpopulated bins are not considered further.
  std::unordered_map<int64_t, std::vector<int>> binLocs;
  for (int gloc = 0; gloc < gnlocs; ++gloc) {
    binLocs[bin[gloc]].push_back(gloc);
  }
  // Total number of populated bins.
  const size_t numPopulatedBins = binLocs.size();

  // Fill modified lat/lon/time and pathsGroupedByLocation.
  // pathsGroupedByLocation will always be of length `gnlocs`,
  // but the modified lat/lon/time may be shorter if there are
  // multiple locations in the same lat-lon bin.
  pathsGroupedByLocation.reserve(gnlocs);
  latsModified.reserve(numPopulatedBins);
  lonsModified.reserve(numPopulatedBins);
  timesModified.reserve(numPopulatedBins);

  // Bin locations that have been encountered in the
  // loop over global locations.
  std::vector<int> binLocsUsed;

  // Loop over global locations.
  for (int gloc = 0; gloc < gnlocs; ++gloc) {
    // Only process each location on its patch rank.
    if (globalPatchRanks[gloc] == rank) {
      // Bin associated with this location.
      const int64_t binLoc = bin[gloc];
      // First location associated with this bin.
      const int firstBinLoc = binLocs[binLoc].front();
      // Determine whether this bin location has been used before.
      const auto binLocValue =
        std::find(binLocsUsed.begin(), binLocsUsed.end(), firstBinLoc);
      if (binLocValue != binLocsUsed.end()) {
        // The bin location has been used before.
        // Extend pathsGroupedByLocation but do not extend
        // modified lat/lon/time.
        const size_t idx = binLocValue - binLocsUsed.begin();
        pathsGroupedByLocation.push_back({idx, idx + 1});
      } else {
        // This bin location has not been used before.
        // Extend pathsGroupedByLocation and modified lat/lon/time.
        latsModified.push_back(globalLatitude[firstBinLoc]);
        lonsModified.push_back(globalLongitude[firstBinLoc]);
        timesModified.push_back(globalTime[firstBinLoc]);
        binLocsUsed.push_back(firstBinLoc);
        const size_t idx = binLocsUsed.size() - 1;
        pathsGroupedByLocation.push_back({idx, idx + 1});
      }
    }
  }

  // Set number of modified lat/lon/time samples.
  updateNumSamples(latsModified.size());
}

}  // namespace ufo
