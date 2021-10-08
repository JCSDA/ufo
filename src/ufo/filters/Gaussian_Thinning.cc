/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/Gaussian_Thinning.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/GaussianThinningParameters.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/EquispacedBinSelectorBase.h"
#include "ufo/utils/GeodesicDistanceCalculator.h"
#include "ufo/utils/MaxNormDistanceCalculator.h"
#include "ufo/utils/metoffice/MetOfficeSort.h"
#include "ufo/utils/RecursiveSplitter.h"
#include "ufo/utils/RoundingEquispacedBinSelector.h"
#include "ufo/utils/SpatialBinSelector.h"
#include "ufo/utils/TruncatingEquispacedBinSelector.h"

namespace ufo {

// -----------------------------------------------------------------------------

Gaussian_Thinning::Gaussian_Thinning(ioda::ObsSpace & obsdb,
                                     const GaussianThinningParameters & params,
                                     std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                     std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, params, flags, obserr), options_(params)
{
  oops::Log::debug() << "Gaussian_Thinning: config = " << options_ << std::endl;
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::applyFilter(const std::vector<bool> & apply,
                                    const Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  ObsAccessor obsAccessor = createObsAccessor();

  std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply, *flags_,
                                     filtervars, options_.thinIfAnyFilterVariablesAreValid.value());

  if (options_.opsCompatibilityMode) {
    // Sort observations by latitude
    const std::vector<float> lat = obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
    metOfficeSort(validObsIds.begin(), validObsIds.end(), [&lat] (size_t id) { return lat[id]; });
  }

  std::vector<float> distancesToBinCenter(validObsIds.size(), 0.f);
  std::unique_ptr<DistanceCalculator> distanceCalculator = makeDistanceCalculator(options_);

  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(
        validObsIds, options_.opsCompatibilityMode);
  groupObservationsByVerticalCoordinate(validObsIds, *distanceCalculator, obsAccessor,
                                        splitter, distancesToBinCenter);
  groupObservationsByTime(validObsIds, *distanceCalculator, obsAccessor,
                          splitter, distancesToBinCenter);
  groupObservationsBySpatialLocation(validObsIds, *distanceCalculator, obsAccessor,
                                     splitter, distancesToBinCenter);

  const std::vector<bool> isThinned = identifyThinnedObservations(
        validObsIds, obsAccessor, splitter, distancesToBinCenter);
  obsAccessor.flagRejectedObservations(isThinned, flagged);
}

// -----------------------------------------------------------------------------

ObsAccessor Gaussian_Thinning::createObsAccessor() const {
  if (options_.categoryVariable.value() != boost::none) {
    return ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
          obsdb_, *options_.categoryVariable.value());
  } else {
    return ObsAccessor::toAllObservations(obsdb_);
  }
}

// -----------------------------------------------------------------------------

std::unique_ptr<DistanceCalculator> Gaussian_Thinning::makeDistanceCalculator(
    const GaussianThinningParameters &options) {
  DistanceNorm distanceNorm = options.distanceNorm.value().value_or(DistanceNorm::GEODESIC);
  if (options.opsCompatibilityMode)
    distanceNorm = DistanceNorm::MAXIMUM;
  switch (distanceNorm) {
  case DistanceNorm::GEODESIC:
    return std::unique_ptr<DistanceCalculator>(new GeodesicDistanceCalculator());
  case DistanceNorm::MAXIMUM:
    return std::unique_ptr<DistanceCalculator>(new MaxNormDistanceCalculator());
  }
  throw eckit::BadParameter("Unrecognized distance norm", Here());
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::groupObservationsBySpatialLocation(
    const std::vector<size_t> &validObsIds,
    const DistanceCalculator &distanceCalculator,
    const ObsAccessor &obsAccessor,
    RecursiveSplitter &splitter,
    std::vector<float> &distancesToBinCenter) const {
  boost::optional<SpatialBinSelector> binSelector = makeSpatialBinSelector(options_);
  if (binSelector == boost::none)
    return;

  oops::Log::debug() << "Gaussian_Thinning: zonal band width (degrees) = "
                     << binSelector->latitudeBinWidth() << std::endl;
  oops::Log::debug() << "Gaussian_Thinning: number of horizontal bins = "
                     << binSelector->totalNumBins() << std::endl;

  std::vector<float> lat = obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
  std::vector<float> lon = obsAccessor.getFloatVariableFromObsSpace("MetaData", "longitude");
  if (!options_.opsCompatibilityMode) {
    // Longitudes will typically be either in the [-180, 180] degree range or in the [0, 360]
    // degree range. When the OPS compatibility mode is off, the spatial bin selector is constructed
    // with the latter convention in mind, so we need to shift any negative longitudes up by 360
    // degrees.
    for (float &longitude : lon)
      if (longitude < 0)
        longitude += 360;
  }

  std::vector<int> latBins;
  std::vector<int> lonBins;
  latBins.reserve(validObsIds.size());
  lonBins.reserve(validObsIds.size());
  for (size_t obsId : validObsIds) {
    const size_t latBin = binSelector->latitudeBin(lat[obsId]);
    latBins.push_back(latBin);
    lonBins.push_back(binSelector->longitudeBin(latBin, lon[obsId]));
  }
  splitter.groupBy(latBins);
  splitter.groupBy(lonBins);

  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    const size_t obsId = validObsIds[validObsIndex];
    float component = distanceCalculator.spatialDistanceComponent(
          lat[obsId], lon[obsId],
          binSelector->latitudeBinCenter(latBins[validObsIndex]),
          binSelector->longitudeBinCenter(latBins[validObsIndex], lonBins[validObsIndex]),
          binSelector->inverseLatitudeBinWidth(),
          binSelector->inverseLongitudeBinWidth(latBins[validObsIndex]));
    distancesToBinCenter[validObsIndex] = distanceCalculator.combineDistanceComponents(
          distancesToBinCenter[validObsIndex], component);
  }
}

// -----------------------------------------------------------------------------

boost::optional<SpatialBinSelector> Gaussian_Thinning::makeSpatialBinSelector(
    const GaussianThinningParameters &options) {
  if (options.horizontalMesh <= 0)
    return boost::none;

  oops::Log::debug() << "Gaussian_Thinning: requested horizontal bin size (km) = "
                     << options.horizontalMesh << std::endl;

  bool roundHorizontalBinCountToNearest =
      options.roundHorizontalBinCountToNearest.value().value_or(false);
  if (options.opsCompatibilityMode)
    roundHorizontalBinCountToNearest = true;
  SpatialBinCountRoundingMode roundingMode = roundHorizontalBinCountToNearest ?
        SpatialBinCountRoundingMode::NEAREST : SpatialBinCountRoundingMode::DOWN;

  const float earthRadius = Constants::mean_earth_rad;  // km
  const float meridianLength = M_PI * earthRadius;
  const float tentativeNumLatBins = meridianLength / options.horizontalMesh;
  const int numLatBins = SpatialBinSelector::roundNumBins(tentativeNumLatBins, roundingMode);

  if (options.useReducedHorizontalGrid) {
    // Use fewer bins at high latitudes
    return SpatialBinSelector(numLatBins, roundingMode, options.opsCompatibilityMode);
  } else {
    // Use the same number of bins at all latitudes
    const int equatorToMeridianLengthRatio = 2;
    return SpatialBinSelector(numLatBins, equatorToMeridianLengthRatio * numLatBins,
                              options.opsCompatibilityMode);
  }
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::groupObservationsByVerticalCoordinate(
    const std::vector<size_t> &validObsIds,
    const DistanceCalculator &distanceCalculator,
    const ObsAccessor &obsAccessor,
    RecursiveSplitter &splitter,
    std::vector<float> &distancesToBinCenter) const {
  std::unique_ptr<EquispacedBinSelectorBase> binSelector = makeVerticalBinSelector(options_);
  if (!binSelector)
    return;

  if (binSelector->numBins() != boost::none)
    oops::Log::debug() << "Gaussian_Thinning: number of vertical bins = "
                       << *binSelector->numBins() << std::endl;

  std::vector<float> vcoord = obsAccessor.getFloatVariableFromObsSpace(
        "MetaData", options_.verticalCoord);

  std::vector<int> bins;
  bins.reserve(validObsIds.size());
  for (size_t obsId : validObsIds)
  {
    bins.push_back(binSelector->bin(vcoord[obsId]));
  }
  splitter.groupBy(bins);

  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    const size_t obsId = validObsIds[validObsIndex];
    const float component = distanceCalculator.nonspatialDistanceComponent(
          vcoord[obsId], binSelector->binCenter(bins[validObsIndex]),
          binSelector->inverseBinWidth());
    distancesToBinCenter[validObsIndex] = distanceCalculator.combineDistanceComponents(
          distancesToBinCenter[validObsIndex], component);
  }
}

// -----------------------------------------------------------------------------

std::unique_ptr<EquispacedBinSelectorBase> Gaussian_Thinning::makeVerticalBinSelector(
    const GaussianThinningParameters &options) {
  if (options.verticalMesh <= 0)
    return nullptr;

  if (options.opsCompatibilityMode) {
    return std::make_unique<RoundingEquispacedBinSelector>(
          options.verticalMesh, options.verticalMin + options.verticalMesh / 2);
  } else {
    const int numVerticalBins = std::max(
          1,
          static_cast<int>(std::ceil((options.verticalMax - options.verticalMin) /
                                     options.verticalMesh)));
    // Adjust verticalMax upwards to make the range of vertical coordinates
    // evenly divisible into bins of width verticalMesh.
    const float adjustedVerticalMax = options.verticalMin + numVerticalBins * options.verticalMesh;

    oops::Log::debug() << "Gaussian_Thinning: number of vertical bins = "
                     << numVerticalBins << std::endl;
    return std::make_unique<TruncatingEquispacedBinSelector>(
          options.verticalMin, adjustedVerticalMax, numVerticalBins);
  }
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::groupObservationsByTime(
    const std::vector<size_t> &validObsIds,
    const DistanceCalculator &distanceCalculator,
    const ObsAccessor &obsAccessor,
    RecursiveSplitter &splitter,
    std::vector<float> &distancesToBinCenter) const {
  util::DateTime timeOffset;
  std::unique_ptr<EquispacedBinSelectorBase> binSelector =
      makeTimeBinSelector(options_, obsdb_.windowStart(), obsdb_.windowEnd(), timeOffset);
  if (!binSelector)
    return;

  if (binSelector->numBins() != boost::none)
    oops::Log::debug() << "Gaussian_Thinning: number of time bins = "
                       << *binSelector->numBins() << std::endl;

  std::vector<util::DateTime> times = obsAccessor.getDateTimeVariableFromObsSpace(
        "MetaData", "datetime");

  std::vector<int> bins;
  bins.reserve(validObsIds.size());
  for (size_t obsId : validObsIds)
  {
    bins.push_back(binSelector->bin((times[obsId] - timeOffset).toSeconds()));
  }
  splitter.groupBy(bins);

  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    const size_t obsId = validObsIds[validObsIndex];
    const float component = distanceCalculator.nonspatialDistanceComponent(
          (times[obsId] - timeOffset).toSeconds(),
          binSelector->binCenter(bins[validObsIndex]),
          binSelector->inverseBinWidth());
    distancesToBinCenter[validObsIndex] = distanceCalculator.combineDistanceComponents(
          distancesToBinCenter[validObsIndex], component);
  }
}

// -----------------------------------------------------------------------------

std::unique_ptr<EquispacedBinSelectorBase> Gaussian_Thinning::makeTimeBinSelector(
    const GaussianThinningParameters &options,
    const util::DateTime &windowStart,
    const util::DateTime &windowEnd,
    util::DateTime &timeOffset) {
  if (options.timeMesh.value() == boost::none ||
      options.timeMin.value() == boost::none ||
      options.timeMax.value() == boost::none)
    return nullptr;

  const util::Duration timeMesh = options.timeMesh.value().get();
  if (timeMesh.toSeconds() == 0)
    return nullptr;

  const util::DateTime timeMin = options.timeMin.value().get();
  const util::DateTime timeMax = options.timeMax.value().get();

  oops::Log::debug() << "(timeMax - timeMin).toSeconds() = "
                     << ((timeMax - timeMin).toSeconds()) << std::endl;
  oops::Log::debug() << "timeMesh.toSeconds() = " << timeMesh.toSeconds() << std::endl;

  timeOffset = timeMin;

  if (options.opsCompatibilityMode) {
    // Put bin 0 at the center of the assimilation window.
    const util::Duration windowLength = windowEnd - windowStart;
    const int numFullBinsLeftOfWindowCenter = (windowLength / 2).toSeconds() / timeMesh.toSeconds();
    const util::DateTime bin0Center = timeMin + numFullBinsLeftOfWindowCenter * timeMesh +
                                      timeMesh / 2;
    return std::make_unique<RoundingEquispacedBinSelector>(
          timeMesh.toSeconds(), (bin0Center - timeOffset).toSeconds());
  } else {
    // The upper bound of the time interval is effectively adjusted upwards
    // to make space for an integral number of bins of specified width.
    const int numTimeBins = std::max(
          1,
          static_cast<int>(std::ceil((timeMax - timeMin).toSeconds() /
                                     static_cast<float>(timeMesh.toSeconds()))));

    return std::make_unique<TruncatingEquispacedBinSelector>(
          0.0f, numTimeBins * timeMesh.toSeconds(), numTimeBins);
  }
}

// -----------------------------------------------------------------------------

std::vector<bool> Gaussian_Thinning::identifyThinnedObservations(
    const std::vector<size_t> &validObsIds,
    const ObsAccessor &obsAccessor,
    const RecursiveSplitter &splitter,
    const std::vector<float> &distancesToBinCenter) const {

  std::function<bool(size_t, size_t)> comparator = makeObservationComparator(
        validObsIds, distancesToBinCenter, obsAccessor);

  size_t totalNumObs = obsAccessor.totalNumObservations();

  std::vector<bool> isThinned(totalNumObs, false);
  for (auto group : splitter.multiElementGroups()) {
    const size_t bestValidObsIndex = *std::min_element(
          std::begin(group), std::end(group), comparator);

    for (size_t validObsIndex : group)
      if (validObsIndex != bestValidObsIndex)
        isThinned[validObsIds[validObsIndex]] = true;
  }

  return isThinned;
}

// Should return true if the first observation is "better" than the second one.
std::function<bool(size_t, size_t)> Gaussian_Thinning::makeObservationComparator(
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &distancesToBinCenter,
    const ObsAccessor &obsAccessor) const
{
  if (options_.priorityVariable.value() == boost::none) {
    return [&distancesToBinCenter](size_t validObsIndexA, size_t validObsIndexB) {
      return distancesToBinCenter[validObsIndexA] < distancesToBinCenter[validObsIndexB];
    };
  }

  const ufo::Variable priorityVariable = options_.priorityVariable.value().get();

  std::vector<int> priorities = obsAccessor.getIntVariableFromObsSpace(
        priorityVariable.group(), priorityVariable.variable());

  // TODO(wsmigaj): In C++14, use move capture for 'priorities'.
  return [priorities, &validObsIds, &distancesToBinCenter]
         (size_t validObsIndexA, size_t validObsIndexB) {
      // Prefer observations with large priorities and small distances
      return std::make_pair(-priorities[validObsIds[validObsIndexA]],
                            distancesToBinCenter[validObsIndexA]) <
             std::make_pair(-priorities[validObsIds[validObsIndexB]],
                            distancesToBinCenter[validObsIndexB]);
    };
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::print(std::ostream & os) const {
  os << "Gaussian_Thinning: config = " << options_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
