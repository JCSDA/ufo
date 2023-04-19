/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficePressureConsistencyCheck.h"

#include <algorithm>

#include "ufo/filters/ObsAccessor.h"
#include "ufo/utils/metoffice/MetOfficeSort.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

MetOfficePressureConsistencyCheck::MetOfficePressureConsistencyCheck
(ioda::ObsSpace & obsdb,
 const Parameters_ &parameters,
 std::shared_ptr<ioda::ObsDataVector<int> > flags,
 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::debug() << "MetOfficePressureConsistencyCheck: config = " << options_ << std::endl;
}

// Required for the correct destruction of options_.
MetOfficePressureConsistencyCheck::~MetOfficePressureConsistencyCheck()
{}

void MetOfficePressureConsistencyCheck::applyFilter(const std::vector<bool> & apply,
                                          const Variables & filtervars,
                                          std::vector<std::vector<bool>> & flagged) const {
    ObsAccessor obsAccessor = createObsAccessor();

    const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply, *flags_, filtervars);

    RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);

    std::vector<util::DateTime> times = obsAccessor.getDateTimeVariableFromObsSpace(
          "MetaData", "dateTime");
    std::vector<bool> PmslUsedFlag = obsAccessor.getBoolVariableFromObsSpace(
        "DiagnosticFlags/PmslUsed", "surfacePressure");
    std::vector<bool> PstdUsedFlag = obsAccessor.getBoolVariableFromObsSpace(
        "DiagnosticFlags/PstdUsed", "surfacePressure");
    std::vector<bool> PstnUsedFlag = obsAccessor.getBoolVariableFromObsSpace(
        "DiagnosticFlags/PstnUsed", "surfacePressure");

    std::vector<bool> isRejected(times.size(), false);

    splitter.sortGroupsBy([&times, &validObsIds](size_t obsIndex)
                          { return times[validObsIds[obsIndex]]; });

    for (RecursiveSplitter::Group group : splitter.multiElementGroups()) {
      // First find the time of the main observation
      std::vector<size_t>::const_iterator seedIt;
      if (options_.seedTime.value() == boost::none) {
        seedIt = group.begin();
      } else {
        seedIt = findSeed(validObsIds,
              times, group.begin(), group.end(), *options_.seedTime.value());
      }

      // Then find the originating pressure source of the main observation
      bool PmslUsed = false;
      bool PstdUsed = false;
      bool PstnUsed = false;

      if (PmslUsedFlag[validObsIds[*seedIt]] == 1) {PmslUsed = true;}
      else if (PstdUsedFlag[validObsIds[*seedIt]] == 1) {PstdUsed = true;}
      else if (PstnUsedFlag[validObsIds[*seedIt]] == 1) {PstnUsed = true;}

      // Then check that all other accepted observations in the group use the same originating
      // pressure souce as the main observation
      for (std::vector<size_t>::const_iterator it = group.begin(); it != group.end(); ++it) {
        const size_t validObsIndex = *it;

        if (PmslUsed) {
          if (not PmslUsedFlag[validObsIds[validObsIndex]])
            isRejected[validObsIds[validObsIndex]] = true;
        } else if (PstnUsed) {
          if (not PstnUsedFlag[validObsIds[validObsIndex]])
            isRejected[validObsIds[validObsIndex]] = true;
        } else if (PstdUsed) {
          if (not PstdUsedFlag[validObsIds[validObsIndex]])
            isRejected[validObsIds[validObsIndex]] = true;
        }
      }
    }

    // Finally reject any of the obsrvations where the originating pressure souce differes
    // from the main observation
    obsAccessor.flagRejectedObservations(isRejected, flagged);
}

std::vector<size_t>::const_iterator MetOfficePressureConsistencyCheck::findSeed(
    std::vector<size_t> validObsIds,
    std::vector<util::DateTime> times,
    std::vector<size_t>::const_iterator validObsIndicesBegin,
    std::vector<size_t>::const_iterator validObsIndicesEnd,
    util::DateTime targetTime) const {
  ASSERT_MSG(validObsIndicesEnd - validObsIndicesBegin != 0,
             "The range of observation indices must not be empty");

  auto isEarlierThan = [&](size_t validObsIndexA, const util::DateTime &timeB) {
    return (times[validObsIds[validObsIndexA]] < timeB);
  };


  const std::vector<size_t>::const_iterator firstGreaterOrEqualToTargetIt =
      std::lower_bound(validObsIndicesBegin, validObsIndicesEnd, targetTime, isEarlierThan);
  if (firstGreaterOrEqualToTargetIt == validObsIndicesBegin) {
    return firstGreaterOrEqualToTargetIt;
  }
  if (firstGreaterOrEqualToTargetIt == validObsIndicesEnd) {
    // All observations were taken before targetTime. Return the iterator pointing to the
    // last valid element of the range.
    return validObsIndicesEnd - 1;
  }

  const std::vector<size_t>::const_iterator lastLessThanTargetIt =
      firstGreaterOrEqualToTargetIt - 1;

  util::DateTime firstGreaterOrEqualToTargetTime = times[validObsIds[*firstGreaterOrEqualToTargetIt]];
  util::DateTime lastLessThanTargetTime = times[validObsIds[*lastLessThanTargetIt]];

  // Prefer the later observation if there's a tie
  if ((firstGreaterOrEqualToTargetTime - targetTime).toSeconds() <=
      (targetTime - lastLessThanTargetTime).toSeconds()) {
    return firstGreaterOrEqualToTargetIt;
  } else {
    return lastLessThanTargetIt;
  }
}

ObsAccessor MetOfficePressureConsistencyCheck::createObsAccessor() const {
  return ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
            obsdb_, *options_.categoryVariable.value() );
}

void MetOfficePressureConsistencyCheck::print(std::ostream & os) const {
  os << "MetOfficePressureConsistencyCheck: config = " << options_ << std::endl;
}

}  // namespace ufo
