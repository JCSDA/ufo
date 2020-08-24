/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/PoissonDiskThinning.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/container/KDTree.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/IsAnyPointInVolumeInterior.h"
#include "oops/util/Logger.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/PoissonDiskThinningParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

namespace {

/// \brief Abstract interface of a container storing point sets and able to answer spatial queries
/// needed by PoissonDiskThinning ("does any point lie in the interior of an axis-aligned
/// ellipsoid/cylinder?").
template <int numDims_>
class PointIndex {
 public:
  typedef float CoordType;
  static const int numDims = numDims_;

  typedef std::array<CoordType, numDims> Point;
  typedef std::array<CoordType, numDims> Extent;

  virtual ~PointIndex() {}

  virtual void insert(const Point &point) = 0;

  virtual bool isAnyPointInCylinderInterior(const Point &center,
                                            const Extent &semiAxes,
                                            int numSpatialDims) const = 0;

  virtual bool isAnyPointInEllipsoidInterior(const Point &center,
                                             const Extent &semiAxes) const = 0;
};

/// \brief An implementation of PointIndex storing the point set in a kd-tree.
template <int numDims_>
class KDTree : public PointIndex<numDims_> {
 public:
  typedef PointIndex<numDims_> Base;

  typedef typename Base::CoordType CoordType;
  typedef typename Base::Point Point;
  typedef typename Base::Extent Extent;

  static const int numDims = Base::numDims;

  void insert(const Point &point) override;

  bool isAnyPointInCylinderInterior(const Point &center,
                                    const Extent &semiAxes,
                                    int numSpatialDims) const override;

  bool isAnyPointInEllipsoidInterior(const Point &center,
                                     const Extent &semiAxes) const override;

 private:
  struct EmptyPayload {};

  struct TreeTraits {
    typedef eckit::geometry::KPoint<numDims> Point;
    typedef EmptyPayload Payload;
  };

  typedef eckit::KDTreeMemory<TreeTraits> KDTreeImpl;
  typedef typename KDTreeImpl::Alloc Alloc;
  typedef typename KDTreeImpl::Node Node;
  typedef typename KDTreeImpl::Point KPoint;
  typedef typename KDTreeImpl::Value Value;

  KDTreeImpl tree_;
};

template <int numDims_>
void KDTree<numDims_>::insert(
    const Point &point) {
  tree_.insert(Value(KPoint(point), EmptyPayload()));
}

template <int numDims_>
bool KDTree<numDims_>::isAnyPointInEllipsoidInterior(
    const Point &center, const Extent &semiAxes) const {
  KPoint lbound, ubound;
  for (int d = 0; d < numDims; ++d) {
    lbound.data()[d] = center[d] - semiAxes[d];
    ubound.data()[d] = center[d] + semiAxes[d];
  }
  return util::isAnyPointInEllipsoidInterior(tree_, lbound, ubound);
}

template <int numDims_>
bool KDTree<numDims_>::isAnyPointInCylinderInterior(
    const Point &center, const Extent &semiAxes, int numSpatialDims) const {
  KPoint lbound, ubound;
  for (int d = 0; d < numDims; ++d) {
    lbound.data()[d] = center[d] - semiAxes[d];
    ubound.data()[d] = center[d] + semiAxes[d];
  }
  return util::isAnyPointInCylinderInterior(tree_, lbound, ubound, numSpatialDims);
}

}  // namespace

struct PoissonDiskThinning::ObsData
{
  boost::optional<util::ScalarOrMap<int, float>> minHorizontalSpacings;
  boost::optional<std::vector<float>> latitudes;
  boost::optional<std::vector<float>> longitudes;

  boost::optional<util::ScalarOrMap<int, float>> minVerticalSpacings;
  boost::optional<std::vector<float>> pressures;

  boost::optional<util::ScalarOrMap<int, util::Duration>> minTimeSpacings;
  boost::optional<std::vector<util::DateTime>> times;

  boost::optional<std::vector<int>> priorities;
};

PoissonDiskThinning::PoissonDiskThinning(ioda::ObsSpace & obsdb,
                                         const eckit::Configuration & config,
                                         std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                         std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "PoissonDiskThinning: config = " << config_ << std::endl;

  options_.reset(new PoissonDiskThinningParameters());
  options_->deserialize(config);
}

// Required for the correct destruction of options_.
PoissonDiskThinning::~PoissonDiskThinning()
{}

void PoissonDiskThinning::applyFilter(const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  const std::vector<size_t> validObsIds = getValidObservationIds(apply);

  if (validObsIds.empty()) {
    return;
  }

  int numSpatialDims, numNonspatialDims;
  ObsData obsData = getObsData(numSpatialDims, numNonspatialDims);

  std::vector<bool> isThinned(apply.size(), false);

  // Thin points from each category separately.
  RecursiveSplitter categorySplitter(validObsIds.size());
  groupObservationsByCategory(validObsIds, categorySplitter);
  for (auto categoryGroup : categorySplitter.multiElementGroups()) {
    std::vector<size_t> obsIdsInCategory;
    for (size_t validObsIndex : categoryGroup) {
      obsIdsInCategory.push_back(validObsIds[validObsIndex]);
    }

    // Within each category, sort points by descending priority and then (if requested)
    // randomly shuffle points of equal priority.
    RecursiveSplitter prioritySplitter(obsIdsInCategory.size());
    groupObservationsByPriority(obsIdsInCategory, prioritySplitter);
    if (options_->shuffle) {
      if (options_->randomSeed.value() != boost::none)
        prioritySplitter.shuffleGroups(*options_->randomSeed.value());
      else
        prioritySplitter.shuffleGroups();
    }

    // Select points to retain within the category.
    thinCategory(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, numNonspatialDims,
                 isThinned);
  }

  flagThinnedObservations(isThinned, flagged);

  if (filtervars.size() != 0) {
    oops::Log::trace() << "PoissonDiskThinning: flagged? = " << flagged[0] << std::endl;
  }
}

PoissonDiskThinning::ObsData PoissonDiskThinning::getObsData(int &numSpatialDims,
                                                             int &numNonspatialDims) const
{
  ObsData obsData;

  numSpatialDims = 0;
  numNonspatialDims = 0;

  {
    obsData.minHorizontalSpacings = options_->minHorizontalSpacing.value();
    if (obsData.minHorizontalSpacings != boost::none) {
      validateSpacings(*obsData.minHorizontalSpacings, "min_horizontal_spacing");
      obsData.latitudes.emplace(obsdb_.nlocs());
      obsData.longitudes.emplace(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "latitude", *obsData.latitudes);
      obsdb_.get_db("MetaData", "longitude", *obsData.longitudes);
      numSpatialDims = 3;
    }
  }

  {
    obsData.minVerticalSpacings = options_->minVerticalSpacing.value();
    if (obsData.minVerticalSpacings != boost::none) {
      validateSpacings(*obsData.minVerticalSpacings, "min_vertical_spacing");
      obsData.pressures.emplace(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "air_pressure", *obsData.pressures);
      ++numNonspatialDims;
    }
  }

  {
    obsData.minTimeSpacings = options_->minTimeSpacing.value();
    if (obsData.minTimeSpacings != boost::none) {
      validateSpacings(*obsData.minTimeSpacings, "min_time_spacing");
      obsData.times.emplace(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "datetime", *obsData.times);
      ++numNonspatialDims;
    }
  }

  {
    const boost::optional<Variable> priorityVariable = options_->priorityVariable;
    if (priorityVariable != boost::none) {
      obsData.priorities.emplace(obsdb_.nlocs());
      obsdb_.get_db(priorityVariable.get().group(), priorityVariable.get().variable(),
                    *obsData.priorities);
    }
  }

  return obsData;
}

template <typename ValueType>
void PoissonDiskThinning::validateSpacings(
    const util::ScalarOrMap<int, ValueType> &spacingsByPriority,
    const std::string &parameterName) const {
  if (spacingsByPriority.isScalar())
    return;

  if (spacingsByPriority.begin() == spacingsByPriority.end())
    throw eckit::BadParameter(parameterName + " must be a scalar or a non-empty map");

  // The map is ordered by increasing priority, so the spacing of every item must be
  // no larger than that of the previous item
  ValueType prevSpacing = spacingsByPriority.begin()->second;
  for (const auto &priorityAndSpacing : spacingsByPriority) {
    const ValueType &spacing = priorityAndSpacing.second;
    if (spacing > prevSpacing)
      throw eckit::BadParameter(parameterName +
                                ": exclusion volumes of lower-priority observations must be "
                                "at least as large as those of higher-priority ones.");
  }
}

std::vector<size_t> PoissonDiskThinning::getValidObservationIds(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && (*flags_)[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsId);
  return validObsIds;
}

void PoissonDiskThinning::groupObservationsByCategory(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  boost::optional<Variable> categoryVariable = options_->categoryVariable;
  if (categoryVariable == boost::none)
    return;

  ioda::ObsDataVector<int> obsDataVector(obsdb_, categoryVariable.get().variable(),
                                         categoryVariable.get().group());
  const ioda::ObsDataRow<int> &category = obsDataVector[0];

  std::vector<int> validObsCategories(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex)
    validObsCategories[validObsIndex] = category[validObsIds[validObsIndex]];
  splitter.groupBy(validObsCategories);
}

void PoissonDiskThinning::groupObservationsByPriority(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  boost::optional<Variable> priorityVariable = options_->priorityVariable;
  if (priorityVariable == boost::none)
    return;

  ioda::ObsDataVector<int> obsDataVector(obsdb_, priorityVariable.get().variable(),
                                         priorityVariable.get().group());
  const ioda::ObsDataRow<int> &priority = obsDataVector[0];

  auto reverse = [](int i) {
      return -i - std::numeric_limits<int>::lowest() + std::numeric_limits<int>::max();
  };

  std::vector<int> validObsPriorities(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex)
    // reversing because we want to start with the highest-priority items
    validObsPriorities[validObsIndex] = reverse(priority[validObsIds[validObsIndex]]);
  splitter.groupBy(validObsPriorities);
}

void PoissonDiskThinning::thinCategory(const ObsData &obsData,
                                       const std::vector<size_t> &obsIdsInCategory,
                                       const RecursiveSplitter &prioritySplitter,
                                       int numSpatialDims,
                                       int numNonspatialDims,
                                       std::vector<bool> &isThinned) const {
  switch (numSpatialDims + numNonspatialDims) {
  case 0:
    return;  // nothing to do
  case 1:
    return thinCategory<1>(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, isThinned);
  case 2:
    return thinCategory<2>(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, isThinned);
  case 3:
    return thinCategory<3>(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, isThinned);
  case 4:
    return thinCategory<4>(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, isThinned);
  case 5:
    return thinCategory<5>(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims, isThinned);
  }

  ABORT("Unexpected number of thinning dimensions");
}

template <int numDims>
void PoissonDiskThinning::thinCategory(const ObsData &obsData,
                                       const std::vector<size_t> &obsIdsInCategory,
                                       const RecursiveSplitter &prioritySplitter,
                                       int numSpatialDims,
                                       std::vector<bool> &isThinned) const {
  KDTree<numDims> pointIndex;

  for (auto priorityGroup : prioritySplitter.groups()) {
    for (size_t obsIndex : priorityGroup) {
      const size_t obsId = obsIdsInCategory[obsIndex];
      std::array<float, numDims> point = getObservationPosition<numDims>(obsId, obsData);
      std::array<float, numDims> semiAxes = getExclusionVolumeSemiAxes<numDims>(obsId, obsData);
      if ((options_->exclusionVolumeShape == ExclusionVolumeShape::CYLINDER &&
           pointIndex.isAnyPointInCylinderInterior(point, semiAxes, numSpatialDims)) ||
          (options_->exclusionVolumeShape == ExclusionVolumeShape::ELLIPSOID &&
           pointIndex.isAnyPointInEllipsoidInterior(point, semiAxes))) {
        isThinned[obsId] = true;
      } else {
        pointIndex.insert(point);
      }
    }
  }
}

template <int numDims>
std::array<float, numDims> PoissonDiskThinning::getObservationPosition(
    size_t obsId, const ObsData &obsData) const {
  std::array<float, numDims> position;

  unsigned int dim = 0;

  if (obsData.latitudes && obsData.longitudes) {
    const float deg2rad = static_cast<float>(M_PI / 180.0);
    const float earthRadius = Constants::mean_earth_rad;

    const float lon = deg2rad * (*obsData.longitudes)[obsId];
    const float lat = deg2rad * (*obsData.latitudes)[obsId];
    const float sinLat = std::sin(lat);
    const float cosLat = std::cos(lat);
    const float sinLon = std::sin(lon);
    const float cosLon = std::cos(lon);

    position[dim++] = earthRadius * cosLat * cosLon;
    position[dim++] = earthRadius * cosLat * sinLon;
    position[dim++] = earthRadius * sinLat;
  }

  if (obsData.pressures) {
    position[dim++] = (*obsData.pressures)[obsId];
  }

  if (obsData.times) {
    // We choose obsData.times[0] as the reference time when converting datetimes to floats.
    // Maybe there's a better way.
    position[dim++] = ((*obsData.times)[obsId] - (*obsData.times)[0]).toSeconds();
  }

  return position;
}

template <int numDims>
std::array<float, numDims> PoissonDiskThinning::getExclusionVolumeSemiAxes(
    size_t obsId, const ObsData &obsData) const {

  std::array<float, numDims> semiAxes;

  const int priority = obsData.priorities == boost::none ? 0 : (*obsData.priorities)[obsId];

  unsigned int dim = 0;

  if (obsData.minHorizontalSpacings) {
    const float earthDiameter = 2 * Constants::mean_earth_rad;
    const float invEarthDiameter = 1 / earthDiameter;

    const float minGeodesicDistance = obsData.minHorizontalSpacings->at(priority);
    const float minEuclideanDistance =
        earthDiameter * std::sin(minGeodesicDistance * invEarthDiameter);

    semiAxes[dim++] = minEuclideanDistance;
    semiAxes[dim++] = minEuclideanDistance;
    semiAxes[dim++] = minEuclideanDistance;
  }

  if (obsData.minVerticalSpacings) {
    semiAxes[dim++] = obsData.minVerticalSpacings->at(priority);
  }

  if (obsData.minTimeSpacings) {
    semiAxes[dim++] = obsData.minTimeSpacings->at(priority).toSeconds();
  }

  return semiAxes;
}

void PoissonDiskThinning::flagThinnedObservations(
    const std::vector<bool> & isThinned,
    std::vector<std::vector<bool>> & flagged) const {
  for (std::vector<bool> & variableFlagged : flagged)
    for (size_t obsId = 0; obsId < isThinned.size(); ++obsId)
      if (isThinned[obsId])
        variableFlagged[obsId] = true;
}

void PoissonDiskThinning::print(std::ostream & os) const {
  os << "PoissonDiskThinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
