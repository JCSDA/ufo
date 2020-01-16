/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/TemporalThinning.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/filters/TemporalThinningParameters.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

namespace {

struct Observation {
  size_t id;
  util::DateTime time;
  int priority = 0;
};

/// \brief Responsible for the selection of observations to retain.
class TemporalThinner {
 public:
  TemporalThinner(const std::vector<size_t> &validObsIds,
                  const std::vector<util::DateTime> &times,
                  const std::vector<int> *priorities,
                  const RecursiveSplitter &splitter,
                  const TemporalThinningParameters &options);

  std::vector<bool> identifyThinnedObservations(size_t totalNumObservations) const;

 private:
  typedef std::vector<size_t>::const_iterator ForwardValidObsIndexIterator;
  typedef std::vector<size_t>::const_reverse_iterator BackwardValidObsIndexIterator;

  Observation getObservation(size_t validObsIndex) const {
    Observation obs;
    obs.id = validObsIds_[validObsIndex];
    obs.time = times_[obs.id];
    if (priorities_)
      obs.priority = (*priorities_)[obs.id];
    return obs;
  }

  util::DateTime getTime(size_t validObsIndex) const {
    return times_[validObsIds_[validObsIndex]];
  }

  int getPriority(size_t validObsIndex) const {
    assert(priorities_ != nullptr);
    return (*priorities_)[validObsIds_[validObsIndex]];
  }

  bool hasPriorities() const {
    return priorities_ != nullptr;
  }

  /// Thin the valid observations with indices from the specified range. The first observation
  /// to retain must be taken at or after \p deadline.
  void thinRangeForwards(ForwardValidObsIndexIterator validIndicesBegin,
                         ForwardValidObsIndexIterator validIndicesEnd,
                         util::DateTime deadline,
                         std::vector<bool> &isThinned) const;

  /// Thin the valid observations with indices from the specified range. The first observation
  /// to retain must be taken at or before \p deadline.
  void thinRangeBackwards(BackwardValidObsIndexIterator validIndicesBegin,
                          BackwardValidObsIndexIterator validIndicesEnd,
                          util::DateTime deadline,
                          std::vector<bool> &isThinned) const;

  template <typename Iterator, typename IsPast, typename IsAtOrPast, typename Advance>
  void thinRange(Iterator validIndicesBegin,
                 Iterator validIndicesEnd,
                 util::DateTime deadline,
                 IsPast isPast,
                 IsAtOrPast isAtOrPast,
                 Advance advance,
                 std::vector<bool> &isThinned) const;

  /// Return an iterator to the first observation to be retained.
  ForwardValidObsIndexIterator findSeed(
      ForwardValidObsIndexIterator validObsIndicesBegin,
      ForwardValidObsIndexIterator validObsIndicesEnd,
      util::DateTime seedTime) const;

  /// Return an iterator to the valid observation taken at a time closest to \p targetTime.
  /// In case of a tie, the later (more recent) observation is selected.
  ForwardValidObsIndexIterator findNearest(
      ForwardValidObsIndexIterator validObsIndicesBegin,
      ForwardValidObsIndexIterator validObsIndicesEnd,
      util::DateTime targetTime) const;

 private:
  const std::vector<size_t> &validObsIds_;
  const std::vector<util::DateTime> &times_;
  const std::vector<int> *priorities_;
  const RecursiveSplitter &splitter_;
  const TemporalThinningParameters &options_;
};

TemporalThinner::TemporalThinner(const std::vector<size_t> &validObsIds,
                                 const std::vector<util::DateTime> &times,
                                 const std::vector<int> *priorities,
                                 const RecursiveSplitter &splitter,
                                 const TemporalThinningParameters &options) :
  validObsIds_(validObsIds),
  times_(times),
  priorities_(priorities),
  splitter_(splitter),
  options_(options)
{}

std::vector<bool> TemporalThinner::identifyThinnedObservations(size_t totalNumObservations) const {
  std::vector<bool> isThinned(totalNumObservations, false);

  for (RecursiveSplitter::Group group : splitter_.multiElementGroups()) {
    if (options_.seedTime.value() == boost::none) {
      util::DateTime deadline = getTime(*group.begin());
      thinRangeForwards(group.begin(), group.end(), deadline, isThinned);
    } else {
      const ForwardValidObsIndexIterator seedIt = findSeed(
            group.begin(), group.end(), *options_.seedTime.value());
      Observation seed = getObservation(*seedIt);
      {
        util::DateTime deadline = seed.time + options_.minSpacing;
        thinRangeForwards(seedIt + 1, group.end(), deadline, isThinned);
      }
      {
        const BackwardValidObsIndexIterator seedRevIt(seedIt + 1);
        const BackwardValidObsIndexIterator revEnd(group.begin());
        util::DateTime deadline = seed.time - options_.minSpacing;
        thinRangeBackwards(seedRevIt, revEnd, seed.time, isThinned);
      }
    }
  }

  return isThinned;
}

void TemporalThinner::thinRangeForwards(ForwardValidObsIndexIterator validIndicesBegin,
                                        ForwardValidObsIndexIterator validIndicesEnd,
                                        util::DateTime deadline,
                                        std::vector<bool> &isThinned) const {
  thinRange(validIndicesBegin, validIndicesEnd, deadline,
            std::greater<util::DateTime>(),
            std::greater_equal<util::DateTime>(),
            // Unfortunately in C++11 std::plus requires both arguments to be of the same type
            [](const util::DateTime &time, const util::Duration &offset) { return time + offset; },
            isThinned);
}

void TemporalThinner::thinRangeBackwards(
    BackwardValidObsIndexIterator validIndicesBegin,
    BackwardValidObsIndexIterator validIndicesEnd,
    util::DateTime deadline,
    std::vector<bool> &isThinned) const {
  thinRange(validIndicesBegin, validIndicesEnd, deadline,
            std::less<util::DateTime>(),
            std::less_equal<util::DateTime>(),
            [](const util::DateTime &time, const util::Duration &offset) { return time - offset; },
            isThinned);
}

template <typename Iterator, typename IsPast, typename IsAtOrPast, typename Advance>
void TemporalThinner::thinRange(Iterator validIndicesBegin,
                                Iterator validIndicesEnd,
                                util::DateTime deadline,
                                IsPast isPast,
                                IsAtOrPast isAtOrPast,
                                Advance advance,
                                std::vector<bool> &isThinned) const {
  boost::optional<Observation> best;
  for (Iterator it = validIndicesBegin; it != validIndicesEnd; ++it) {
    const size_t validObsIndex = *it;
    Observation current = getObservation(validObsIndex);
    if (best != boost::none) {
      // We're looking for a higher-priority observation at or before the deadline
      if (isPast(current.time, deadline)) {
        // We haven't found one
        deadline = advance(best->time, options_.minSpacing);
        best = boost::none;
        // The decision whether to thin 'current' will be taken in the next if statement
      } else {
        if (current.priority > best->priority) {
          isThinned[best->id] = true;
          best = current;
        } else {
          isThinned[current.id] = true;
        }
      }
    }

    if (best == boost::none) {
      // We're looking for an observation at or after the deadline
      if (isAtOrPast(current.time, deadline)) {
        best = current;
        deadline = advance(best->time, options_.tolerance);
      } else {
        isThinned[current.id] = true;
      }
    }
  }
}

typename TemporalThinner::ForwardValidObsIndexIterator TemporalThinner::findSeed(
    ForwardValidObsIndexIterator validObsIndicesBegin,
    ForwardValidObsIndexIterator validObsIndicesEnd,
    util::DateTime seedTime) const {
  const ForwardValidObsIndexIterator nearestToSeedIt = findNearest(
        validObsIndicesBegin, validObsIndicesEnd, seedTime);
  if (!hasPriorities()) {
    return nearestToSeedIt;
  }

  util::DateTime nearestToSeedTime = getObservation(*nearestToSeedIt).time;

  ForwardValidObsIndexIterator acceptableBegin = std::lower_bound(
        validObsIndicesBegin, nearestToSeedIt,
        nearestToSeedTime - options_.tolerance,
        [&](size_t validObsIndexA, const util::DateTime &timeB)
        { return getTime(validObsIndexA) < timeB; });
  ForwardValidObsIndexIterator acceptableEnd = std::upper_bound(
        acceptableBegin, validObsIndicesEnd,
        nearestToSeedTime + options_.tolerance,
        [&](const util::DateTime &timeA, size_t validObsIndexB)
        { return timeA < getTime(validObsIndexB); });

  // Find the element with highest priority in the acceptable range.
  return std::max_element(acceptableBegin, acceptableEnd,
                          [&](size_t validObsIndexA, size_t validObsIndexB)
                          { return getPriority(validObsIndexA) < getPriority(validObsIndexB); });
}

typename TemporalThinner::ForwardValidObsIndexIterator TemporalThinner::findNearest(
    ForwardValidObsIndexIterator validObsIndicesBegin,
    ForwardValidObsIndexIterator validObsIndicesEnd,
    util::DateTime targetTime) const {
  ASSERT_MSG(validObsIndicesEnd - validObsIndicesBegin != 0,
             "The range of observation indices must not be empty");

  auto isEarlierThan = [&](size_t validObsIndexA, const util::DateTime &timeB) {
    return getTime(validObsIndexA) < timeB;
  };
  const ForwardValidObsIndexIterator firstGreaterOrEqualToTargetIt =
      std::lower_bound(validObsIndicesBegin, validObsIndicesEnd, targetTime, isEarlierThan);
  if (firstGreaterOrEqualToTargetIt == validObsIndicesBegin) {
    return firstGreaterOrEqualToTargetIt;
  }
  if (firstGreaterOrEqualToTargetIt == validObsIndicesEnd) {
    // All observations were taken before targetTime. Return the iterator pointing to the
    // last valid element of the range.
    return validObsIndicesEnd - 1;
  }

  const ForwardValidObsIndexIterator lastLessThanTargetIt =
      firstGreaterOrEqualToTargetIt - 1;

  Observation firstGreaterOrEqualToTarget = getObservation(*firstGreaterOrEqualToTargetIt);
  Observation lastLessThanTarget = getObservation(*lastLessThanTargetIt);

  // Prefer the later observation if there's a tie
  if ((firstGreaterOrEqualToTarget.time - targetTime).toSeconds() <=
      (targetTime - lastLessThanTarget.time).toSeconds()) {
    return firstGreaterOrEqualToTargetIt;
  } else {
    return lastLessThanTargetIt;
  }
}

}  // namespace

TemporalThinning::TemporalThinning(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                   boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                   boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "TemporalThinning: config = " << config_ << std::endl;

  options_.reset(new TemporalThinningParameters());
  options_->deserialize(config);
}

// Required for the correct destruction of options_.
TemporalThinning::~TemporalThinning()
{}

void TemporalThinning::applyFilter(const std::vector<bool> & apply,
                                   const Variables & filtervars,
                                   std::vector<std::vector<bool>> & flagged) const {
  const std::vector<bool> isThinned = identifyThinnedObservations(apply);

  flagThinnedObservations(isThinned, flagged);

  if (filtervars.size() != 0) {
    oops::Log::trace() << "TemporalThinning: flagged? = " << flagged[0] << std::endl;
  }
}

std::vector<bool> TemporalThinning::identifyThinnedObservations(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds = getValidObservationIds(apply);

  RecursiveSplitter splitter(validObsIds.size());
  groupObservationsByCategory(validObsIds, splitter);

  std::vector<util::DateTime> times(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "datetime", times);
  splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndexA, size_t obsIndexB)
                        { return times[validObsIds[obsIndexA]] < times[validObsIds[obsIndexB]]; });

  std::unique_ptr<ioda::ObsDataVector<int>> prioritiesDataVector = getObservationPriorities();
  const ioda::ObsDataRow<int> *priorities =
      prioritiesDataVector ? &(*prioritiesDataVector)[0] : nullptr;

  TemporalThinner thinner(validObsIds, times, priorities, splitter, *options_);
  return thinner.identifyThinnedObservations(apply.size());
}

std::vector<size_t> TemporalThinning::getValidObservationIds(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && (*flags_)[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsId);
  return validObsIds;
}

void TemporalThinning::groupObservationsByCategory(const std::vector<size_t> &validObsIds,
                                                   RecursiveSplitter &splitter) const {
  boost::optional<Variable> categoryVariable = options_->categoryVariable;
  if (categoryVariable == boost::none)
    return;

  ioda::ObsDataVector<int> obsDataVector(obsdb_, categoryVariable.get().variable(),
                                         categoryVariable.get().group());
  ioda::ObsDataRow<int> &category = obsDataVector[0];

  std::vector<int> validObsCategories(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex)
    validObsCategories[validObsIndex] = category[validObsIds[validObsIndex]];
  splitter.groupBy(validObsCategories);
}

std::unique_ptr<ioda::ObsDataVector<int>> TemporalThinning::getObservationPriorities() const {
  std::unique_ptr<ioda::ObsDataVector<int>> priorities;
  if (options_->priorityVariable.value() != boost::none) {
    const ufo::Variable priorityVariable = options_->priorityVariable.value().get();
    priorities.reset(new ioda::ObsDataVector<int>(
                       obsdb_, priorityVariable.variable(), priorityVariable.group()));
  }
  return priorities;
}

void TemporalThinning::flagThinnedObservations(
    const std::vector<bool> & isThinned,
    std::vector<std::vector<bool>> & flagged) const {
  for (std::vector<bool> & variableFlagged : flagged)
    for (size_t obsId = 0; obsId < isThinned.size(); ++obsId)
       if (isThinned[obsId])
        variableFlagged[obsId] = true;
}

void TemporalThinning::print(std::ostream & os) const {
  os << "TemporalThinning: config = " << config_ << std::endl;
}

}  // namespace ufo
