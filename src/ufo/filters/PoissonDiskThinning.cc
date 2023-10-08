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

#include "eckit/container/KDTree.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/IsAnyPointInVolumeInterior.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/ObsAccessor.h"
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

  virtual bool isAnyPointInBoxInterior(const Point &center,
                                       const Extent &halfSideLength) const = 0;
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

  bool isAnyPointInBoxInterior(const Point &center,
                               const Extent &halfSideLength) const override;

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

template <int numDims_>
bool KDTree<numDims_>::isAnyPointInBoxInterior(
    const Point &center, const Extent &halfSideLength) const {
  KPoint lbound, ubound;
  for (int d = 0; d < numDims; ++d) {
    lbound.data()[d] = center[d] - halfSideLength[d];
    ubound.data()[d] = center[d] + halfSideLength[d];
  }
  return util::isAnyPointInBoxInterior(tree_, lbound, ubound);
}

}  // namespace

struct PoissonDiskThinning::ObsData
{
  boost::optional<util::ScalarOrMap<int, float>> minHorizontalSpacings;
  boost::optional<util::ScalarOrMap<int, float>> minLatitudeSpacings;
  boost::optional<util::ScalarOrMap<int, float>> minLongitudeSpacings;
  boost::optional<std::vector<float>> latitudes;
  boost::optional<std::vector<float>> longitudes;

  boost::optional<util::ScalarOrMap<int, float>> minVerticalSpacings;
  boost::optional<std::vector<float>> pressures;

  boost::optional<util::ScalarOrMap<int, util::Duration>> minTimeSpacings;
  boost::optional<std::vector<util::DateTime>> times;

  boost::optional<std::vector<int>> priorities;

  boost::optional<std::vector<float>> obsForMedian;

  // Total number of observations held by all MPI tasks
  size_t totalNumObs = 0;
};

PoissonDiskThinning::PoissonDiskThinning(ioda::ObsSpace & obsdb,
                                         const Parameters_ &parameters,
                                         std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                         std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::debug() << "PoissonDiskThinning: config = " << options_ << std::endl;
  if (options_.sortVertical.value() != boost::none &&
      options_.minVerticalSpacing.value() == boost::none) {
    throw eckit::UserError(
      ": 'sort_vertical' can only be specified if 'min_vertical_spacing' is specified.", Here());
  } else if (options_.sortVertical.value() != boost::none &&
            (options_.sortVertical.value().value() != "ascending" &&
             options_.sortVertical.value().value() != "descending")) {
    throw eckit::UserError(
      ": 'sort_vertical' can only take on string values 'ascending' or 'descending' "
      "(with respect to the pressure coordinate values).", Here());
  } else if (options_.minHorizontalSpacing.value() != boost::none &&
             (options_.minLatitudeSpacing.value() != boost::none ||
              options_.minLongitudeSpacing.value() != boost::none)) {
    throw eckit::UserError(
       ": can only use either minHorizontalSpacing or both minLatitudeSpacing and "
       "minLongitudeSpacing.", Here());
  } else if ((options_.minLatitudeSpacing.value() != boost::none ||
              options_.minLongitudeSpacing.value() != boost::none) &&
             options_.exclusionVolumeShape != ExclusionVolumeShape::BOX) {
    throw eckit::UserError(
          ": can only use exclusion volume shape BOX with minLatitudeSpacing and "
          "minLongitudeSpacing.", Here());
  } else if ((options_.minLatitudeSpacing.value() == boost::none) !=
             (options_.minLongitudeSpacing.value() == boost::none)) {
    throw eckit::UserError(
          ": must use both minLatitudeSpacing and "
          "minLongitudeSpacing.", Here());
  } else if (!options_.selectMedian.value() && options_.writeMedian.value()) {
    throw eckit::UserError(
          ": write median has no effect if select median is not set.", Here());
  } else if (!options_.selectMedian.value() && options_.opsCompatibilityMode.value()) {
    throw eckit::UserError(
          ": ops compatibility mode has no effect if select median is not set.", Here());
  }
}

// Required for the correct destruction of options_.
PoissonDiskThinning::~PoissonDiskThinning()
{}

void PoissonDiskThinning::applyFilter(const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  ObsAccessor obsAccessor = createObsAccessor();

  const std::vector<size_t> validObsIds = getValidObservationIds(apply, filtervars, obsAccessor);

  int numSpatialDims, numNonspatialDims;
  ObsData obsData = getObsData(obsAccessor, numSpatialDims, numNonspatialDims);

  std::vector<bool> isThinned(obsData.totalNumObs, false);

  if (options_.shuffle) {
    // If the same observations will be processed on multiple ranks, they must use random number
    // generators seeded with the same value because they need to reject the same observations. If
    // all ranks process disjoint sets of observations, synchronisation of random number generators
    // isn't necessary (and in fact somewhat undesirable). The sets of observations processed by
    // different ranks will be disjoint if both of the following conditions are met:
    // * the category_variable option is set to the variable used previously to split the
    //   observation space into records
    // * each record is held on a single rank only (i.e. the observation space isn't using the
    //   InefficientDistribution).
    // Unfortunately the second condition can't be checked without breaking the encapsulation
    // of ioda::Distribution (and reintroducing the isDistributed() method).
    //
    // So to stay on the safe side we synchronise the random number generators whenever the shuffle
    // option is selected.
    synchroniseRandomNumberGenerators(obsdb_.comm());
  }

  // Set up for selecting median observation if required.
  std::vector<float> obsForMedian;
  if (options_.selectMedian) {
    obsForMedian = *obsData.obsForMedian;
  }

  // Thin points from each category separately.
  RecursiveSplitter categorySplitter =
      obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  for (auto categoryGroup : categorySplitter.multiElementGroups()) {
    std::vector<size_t> obsIdsInCategory;
    for (size_t validObsIndex : categoryGroup) {
      obsIdsInCategory.push_back(validObsIds[validObsIndex]);
    }

    // Within each category, sort points by descending priority and then (if requested)
    // randomly shuffle points of equal priority. Otherwise, if 'min_vertical_spacing'
    // and 'sort_vertical' are specified, sort according to pressure coordinate.
    RecursiveSplitter prioritySplitter(obsIdsInCategory.size());
    groupObservationsByPriority(obsIdsInCategory, obsAccessor, prioritySplitter);
    if (options_.shuffle) {
      prioritySplitter.shuffleGroups();
    } else if (options_.sortVertical.value() != boost::none &&
               obsData.minVerticalSpacings != boost::none) {
      const std::vector<float> pressures = *obsData.pressures;
      if (options_.sortVertical.value().value() == "ascending") {
        prioritySplitter.sortGroupsBy([&pressures, &obsIdsInCategory](size_t ind)
                                      {return pressures[obsIdsInCategory[ind]];});
      } else if (options_.sortVertical.value().value() == "descending") {
        prioritySplitter.sortGroupsBy([&pressures, &obsIdsInCategory](size_t ind)
                                      {return -1.0*pressures[obsIdsInCategory[ind]];});
      }
    }
    // Select points to retain within the category.
    if (options_.selectMedian.value()) {
      thinCategoryMedian(obsData, obsIdsInCategory, prioritySplitter, numSpatialDims,
                         numNonspatialDims, obsForMedian, isThinned, options_.opsCompatibilityMode);
    } else {
      thinCategory(obsData, obsIdsInCategory, prioritySplitter,
                   numSpatialDims, numNonspatialDims, isThinned);
    }
  }

  obsAccessor.flagRejectedObservations(isThinned, flagged);

  if (options_.writeMedian) {
    std::vector<float> localObs;
    obsdb_.get_db("ObsValue", filtervars_.variable(0).variable(), localObs);
    for (size_t localObsId = 0; localObsId < obsdb_.nlocs(); localObsId++) {
      if (apply[localObsId]) {
        const size_t globalObsId =
          obsdb_.distribution()->globalUniqueConsecutiveLocationIndex(localObsId);
        localObs[localObsId] = obsForMedian[globalObsId];
      }
    }
    obsdb_.put_db("DerivedObsValue", filtervars_.variable(0).variable(), localObs);
  }
}

ObsAccessor PoissonDiskThinning::createObsAccessor() const {
  if (options_.categoryVariable.value() != boost::none) {
    return ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
          obsdb_, *options_.categoryVariable.value());
  } else {
    return ObsAccessor::toAllObservations(obsdb_);
  }
}

PoissonDiskThinning::ObsData PoissonDiskThinning::getObsData(
    const ObsAccessor &obsAccessor, int &numSpatialDims, int &numNonspatialDims) const
{
  ObsData obsData;

  numSpatialDims = 0;
  numNonspatialDims = 0;

  obsData.minHorizontalSpacings = options_.minHorizontalSpacing.value();
  obsData.minLatitudeSpacings = options_.minLatitudeSpacing.value();
  obsData.minLongitudeSpacings = options_.minLongitudeSpacing.value();
  if (obsData.minHorizontalSpacings != boost::none) {
    validateSpacings(*obsData.minHorizontalSpacings, "min_horizontal_spacing");
    obsData.latitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
    obsData.longitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "longitude");
    numSpatialDims = 3;
  } else if (obsData.minLatitudeSpacings != boost::none ||
             obsData.minLongitudeSpacings != boost::none) {
    validateSpacings(*obsData.minLatitudeSpacings, "min_latitude_spacing");
    validateSpacings(*obsData.minLongitudeSpacings, "min_longitude_spacing");
    obsData.latitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "latitude");
    obsData.longitudes = obsAccessor.getFloatVariableFromObsSpace("MetaData", "longitude");
    numSpatialDims = 2;
  }

  obsData.minVerticalSpacings = options_.minVerticalSpacing.value();
  if (obsData.minVerticalSpacings != boost::none) {
    validateSpacings(*obsData.minVerticalSpacings, "min_vertical_spacing");
    obsData.pressures =
      obsAccessor.getFloatVariableFromObsSpace(options_.pressureGroup, options_.pressureCoord);
    ++numNonspatialDims;
  }

  obsData.minTimeSpacings = options_.minTimeSpacing.value();
  if (obsData.minTimeSpacings != boost::none) {
    validateSpacings(*obsData.minTimeSpacings, "min_time_spacing");
    obsData.times = obsAccessor.getDateTimeVariableFromObsSpace("MetaData", "dateTime");
    ++numNonspatialDims;
  }

  const boost::optional<Variable> priorityVariable = options_.priorityVariable;
  if (priorityVariable != boost::none) {
    obsData.priorities = obsAccessor.getIntVariableFromObsSpace(
          priorityVariable.get().group(), priorityVariable.get().variable());
  }

  if (options_.selectMedian) {
    if (filtervars_.size() != 1) {
      throw eckit::UserError("PoissonDiskThinning error: The select median option will only work "
                             "when the filter is used with a single filter variable.", Here());
    }
    obsData.obsForMedian = obsAccessor.getFloatVariableFromObsSpace("ObsValue",
                                                                filtervars_.variable(0).variable());
  }

  obsData.totalNumObs = obsAccessor.totalNumObservations();

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
    const std::vector<bool> & apply,
    const Variables & filtervars,
    const ObsAccessor &obsAccessor) const {
  std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply, *flags_, filtervars);

  if (!options_.shuffle) {
    // The user wants to process observations in fixed (non-random) order. Ensure the filter
    // produces the same results regardless of the number of MPI ranks by ordering the observations
    // to be processed as if we were running in serial: by record ID.
    const std::vector<size_t> recordIds = obsAccessor.getRecordIds();
    std::stable_sort(validObsIds.begin(), validObsIds.end(),
                     [&recordIds](size_t obsIdA, size_t obsIdB)
                     { return recordIds[obsIdA] < recordIds[obsIdB]; });
  }

  return validObsIds;
}

void PoissonDiskThinning::synchroniseRandomNumberGenerators(const eckit::mpi::Comm &comm) const
{
  const size_t rootRank = 0;

  size_t seed;
  if (options_.randomSeed.value() != boost::none) {
    seed = *options_.randomSeed.value();
  } else {
    if (comm.rank() == rootRank)
      // Perhaps oops could provide a function returning this default seed.
      seed = static_cast<std::uint32_t>(std::time(nullptr));
    comm.broadcast(seed, rootRank);
  }

  RecursiveSplitter splitter(1);
  splitter.setSeed(seed, true /* force? */);
}

void PoissonDiskThinning::groupObservationsByPriority(
    const std::vector<size_t> &validObsIds,
    const ObsAccessor &obsAccessor,
    RecursiveSplitter &splitter) const {
  boost::optional<Variable> priorityVariable = options_.priorityVariable;
  if (priorityVariable == boost::none)
    return;

  // TODO(wsmigaj): reuse the priority vector from obsData.
  std::vector<int> priority = obsAccessor.getIntVariableFromObsSpace(
        priorityVariable.get().group(), priorityVariable.get().variable());

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
    return thinCategory<1>(obsData, obsIdsInCategory, prioritySplitter,
                           numSpatialDims, isThinned);
  case 2:
    return thinCategory<2>(obsData, obsIdsInCategory, prioritySplitter,
                           numSpatialDims, isThinned);
  case 3:
    return thinCategory<3>(obsData, obsIdsInCategory, prioritySplitter,
                           numSpatialDims, isThinned);
  case 4:
    return thinCategory<4>(obsData, obsIdsInCategory, prioritySplitter,
                           numSpatialDims, isThinned);
  case 5:
    return thinCategory<5>(obsData, obsIdsInCategory, prioritySplitter,
                           numSpatialDims, isThinned);
  }

  ABORT("Unexpected number of thinning dimensions");
}

void PoissonDiskThinning::thinCategoryMedian(const ObsData &obsData,
                                             const std::vector<size_t> &obsIdsInCategory,
                                             const RecursiveSplitter &prioritySplitter,
                                             int numSpatialDims,
                                             int numNonspatialDims,
                                             std::vector<float> &obsForMedian,
                                             std::vector<bool> &isThinned,
                                             bool opsCompatibilityMode) const {
  switch (numSpatialDims + numNonspatialDims) {
  case 0:
    return;  // nothing to do
  case 1:
    return thinCategoryMedian<1>(obsData, obsIdsInCategory, prioritySplitter,
                                 numSpatialDims, obsForMedian, isThinned, opsCompatibilityMode);
  case 2:
    return thinCategoryMedian<2>(obsData, obsIdsInCategory, prioritySplitter,
                                 numSpatialDims, obsForMedian, isThinned, opsCompatibilityMode);
  case 3:
    return thinCategoryMedian<3>(obsData, obsIdsInCategory, prioritySplitter,
                                 numSpatialDims, obsForMedian, isThinned, opsCompatibilityMode);
  case 4:
    return thinCategoryMedian<4>(obsData, obsIdsInCategory, prioritySplitter,
                                 numSpatialDims, obsForMedian, isThinned, opsCompatibilityMode);
  case 5:
    return thinCategoryMedian<5>(obsData, obsIdsInCategory, prioritySplitter,
                                 numSpatialDims, obsForMedian, isThinned, opsCompatibilityMode);
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
      if ((options_.exclusionVolumeShape == ExclusionVolumeShape::CYLINDER &&
           pointIndex.isAnyPointInCylinderInterior(point, semiAxes, numSpatialDims)) ||
          (options_.exclusionVolumeShape == ExclusionVolumeShape::ELLIPSOID &&
           pointIndex.isAnyPointInEllipsoidInterior(point, semiAxes)) ||
          (options_.exclusionVolumeShape == ExclusionVolumeShape::BOX &&
           pointIndex.isAnyPointInBoxInterior(point, semiAxes))) {
        isThinned[obsId] = true;
      } else {
        pointIndex.insert(point);
      }
    }
  }
}

template <int numDims>
void PoissonDiskThinning::thinCategoryMedian(const ObsData &obsData,
                                             const std::vector<size_t> &obsIdsInCategory,
                                             const RecursiveSplitter &prioritySplitter,
                                             int numSpatialDims,
                                             std::vector<float> &obsForMedian,
                                             std::vector<bool> &isThinned,
                                             bool opsCompatibilityMode) const {
  // Find all observations enclosed in the exclusion volume that have not previously
  // been found and retain the one closest to the median and thin the rest.
  //
  // Please note that the overlap of exclusion volumes of retained observations may
  // mean that the shape and size of the exclusion volume from which the median is
  // calculated is not consistent with what might be expected. It is expected that
  // this code will be deprecated in the future.
  //
  std::vector<bool> isUsed(obsData.totalNumObs, false);
  const float medianMissingValue = util::missingValue<float>();

  // Function to check if a candidate observation is the median (depends on opsCompatibilityMode).
  auto isMedian = [opsCompatibilityMode](float medianObCandidate,
                                         float obsValueMedianBelow,
                                         float obsValueMedianAbove) {
    const bool obsValueMedianBelowMatch = (medianObCandidate == obsValueMedianBelow);
    if (opsCompatibilityMode) return obsValueMedianBelowMatch;
    const bool obsValueMedianAboveMatch = (medianObCandidate == obsValueMedianAbove);
    return obsValueMedianBelowMatch || obsValueMedianAboveMatch;
  };

  for (auto priorityGroup : prioritySplitter.groups()) {
    // Generate vectors of quantities to avoid having to recalculate them
    // repeatedly.
    std::vector<std::array<float, numDims>> allPoints;
    std::vector<std::array<float, numDims>> allSemiAxes;
    std::vector<size_t> allObsIds;
    for (size_t obsIndex : priorityGroup) {
      const size_t obsId = obsIdsInCategory[obsIndex];
      std::array<float, numDims> point = getObservationPosition<numDims>(obsId, obsData);
      std::array<float, numDims> semiAxes = getExclusionVolumeSemiAxes<numDims>(obsId, obsData);
      allPoints.emplace_back(point);
      allSemiAxes.emplace_back(semiAxes);
      allObsIds.emplace_back(obsId);
    }
    // Split the observations into latitude bins to avoid having to
    // check all observations repeatedly.
    std::vector<int> latBins;
    // Set bin size to be approximately the horizontal spacing in degrees with
    // 10 per cent extra in case of observation spacing falling on the boundary.
    // If horizontal spacing is not specified, by default the bin size is set to
    // be a large number which results in all observations being put into either
    // one or two bins.
    float latBinSize = 180.0;
    if (obsData.minHorizontalSpacings != boost::none) {
      // Convert from km to degrees using the approximate length (in km) of one
      // degree at the equator.
      latBinSize = 1.1 * allSemiAxes[0][0] / 111.0;
    } else if (obsData.minLatitudeSpacings != boost::none) {
      latBinSize = 1.1 * allSemiAxes[0][0];
    }
    for (size_t obsIndex : priorityGroup) {
      const size_t obsId = obsIdsInCategory[obsIndex];
      if (obsData.latitudes != boost::none) {
        latBins.emplace_back(static_cast<int>(std::floor((*obsData.latitudes)[obsId]/latBinSize)));
      } else {
        latBins.emplace_back(0);
      }
    }
    const int minLatBins = *std::min_element(latBins.begin(), latBins.end());
    const int maxLatBins = *std::max_element(latBins.begin(), latBins.end());
    const size_t nLatBins = maxLatBins - minLatBins + 1;
    std::vector<size_t> obsBins[nLatBins];
    for (size_t iObsId = 0; iObsId < allObsIds.size(); iObsId++) {
      obsBins[latBins[iObsId]-minLatBins].emplace_back(iObsId);
    }
    // Loop through the observations until one is found that has not previously
    // been used to calculate a median. Other observations within the exclusion
    // volume of that observation that were not previously used are then used
    // to find the median observation.
    for (size_t iObsId = 0; iObsId < allObsIds.size(); iObsId++) {
      const size_t obsId = allObsIds[iObsId];
      if (isUsed[obsId]) {
        continue;
      }
      // Store the information about this observation.
      std::array<float, numDims> point = allPoints[iObsId];
      KDTree<numDims> pointIndex;
      pointIndex.insert(point);
      std::vector<size_t> medianIndices;
      std::vector<float> medianObs;
      medianIndices.emplace_back(obsId);
      medianObs.emplace_back(obsForMedian[obsId]);
      isUsed[obsId] = true;
      size_t nMedianObs = 1;
      // Find additional observations within the exclusion volume of the
      // observation. Only check the latitude bins surrounding the
      // observation to speed up the processing.
      const int latBin = latBins[iObsId] - minLatBins;
      std::vector<size_t> obsToCheck;
      for (int latOffset = -1; latOffset <= 1; latOffset++) {
        const int latBinOffset = latBin + latOffset;
        if ((latBinOffset < 0) ||
            (latBinOffset >= nLatBins)) {
          continue;
        }
        for (size_t jObsId : obsBins[latBinOffset]) {
          obsToCheck.emplace_back(jObsId);
        }
      }
      std::stable_sort(obsToCheck.begin(), obsToCheck.end());
      for (size_t jObsId : obsToCheck) {
        const size_t obsId2 = allObsIds[jObsId];
        if (isUsed[obsId2]) {
          continue;
        }
        std::array<float, numDims> jPoint = allPoints[jObsId];
        std::array<float, numDims> semiAxes = allSemiAxes[jObsId];
        bool isMatch = false;
        if (options_.exclusionVolumeShape == ExclusionVolumeShape::CYLINDER) {
          isMatch = pointIndex.isAnyPointInCylinderInterior(jPoint, semiAxes, numSpatialDims);
        } else if (options_.exclusionVolumeShape == ExclusionVolumeShape::ELLIPSOID) {
          isMatch = pointIndex.isAnyPointInEllipsoidInterior(jPoint, semiAxes);
        } else if (options_.exclusionVolumeShape == ExclusionVolumeShape::BOX) {
          isMatch = pointIndex.isAnyPointInBoxInterior(jPoint, semiAxes);
        }
        if (isMatch) {
          medianIndices.emplace_back(obsId2);
          medianObs.emplace_back(obsForMedian[obsId2]);
          isUsed[obsId2] = true;
          nMedianObs++;
        }
      }
      // If there is more than one observation in the exclusion volume, find
      // the median observation. If there is an even number of observations,
      // keep the first of the central pair. Thin the rest.
      if (nMedianObs > 1) {
        std::vector<float> medianObsSorted(medianObs);
        std::stable_sort(medianObsSorted.begin(), medianObsSorted.end());
        const size_t obsIndexMedianBelow = static_cast<size_t>(std::floor(0.5 * (nMedianObs - 1)));
        const size_t obsIndexMedianAbove = static_cast<size_t>(std::ceil(0.5 * (nMedianObs - 1)));
        const float obsValueMedianBelow = medianObsSorted[obsIndexMedianBelow];
        const float obsValueMedianAbove = medianObsSorted[obsIndexMedianAbove];
        bool foundMedian = false;
        for (size_t medianIndex = 0; medianIndex < nMedianObs; medianIndex++) {
          if (!foundMedian &&
              isMedian(medianObs[medianIndex], obsValueMedianBelow, obsValueMedianAbove)) {
            foundMedian = true;
            obsForMedian[medianIndices[medianIndex]] =
                                                  0.5 * (obsValueMedianBelow + obsValueMedianAbove);
          } else {
            isThinned[medianIndices[medianIndex]] = true;
            obsForMedian[medianIndices[medianIndex]] = medianMissingValue;
          }
        }
      }
    }  // Loop over observations in the priority group.
  }  // Loop over priority groups.
}

template <int numDims>
std::array<float, numDims> PoissonDiskThinning::getObservationPosition(
    size_t obsId, const ObsData &obsData) const {
  std::array<float, numDims> position;

  unsigned int dim = 0;

  if (obsData.latitudes && obsData.longitudes) {
    if (obsData.minLatitudeSpacings && obsData.minLongitudeSpacings) {
      position[dim++] = (*obsData.latitudes)[obsId];
      position[dim++] = (*obsData.longitudes)[obsId];
    } else {
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
  } else if (obsData.minLatitudeSpacings && obsData.minLongitudeSpacings) {
    semiAxes[dim++] = obsData.minLatitudeSpacings->at(priority);
    semiAxes[dim++] = obsData.minLongitudeSpacings->at(priority);
  }

  if (obsData.minVerticalSpacings) {
    semiAxes[dim++] = obsData.minVerticalSpacings->at(priority);
  }

  if (obsData.minTimeSpacings) {
    semiAxes[dim++] = obsData.minTimeSpacings->at(priority).toSeconds();
  }

  return semiAxes;
}

void PoissonDiskThinning::print(std::ostream & os) const {
  os << "PoissonDiskThinning: config = " << options_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
