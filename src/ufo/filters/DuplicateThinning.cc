/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/DuplicateThinning.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/filters/ObsAccessor.h"
#include "ufo/utils/RecordHandler.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

// -----------------------------------------------------------------------------

DuplicateThinning::DuplicateThinning(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                     ioda::ObsDataVector<int> & flags,
                                     ioda::ObsDataVector<float> & obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "DuplicateThinning constructor start" << std::endl;
  if (parameters_.variableNames.value().size() > 0) {
    for (const Variable & var : parameters_.variableNames.value()) {
      if (obsdb_.has(var.group(), var.variable())) {
        testVars_ +=  var;
      } else {
        const std::string error = "The variable " + var.fullName() + " does not exist.";
        throw eckit::UserError(error, Here());
      }
    }
  }
  if (parameters_.tolerance.value().size() > 0) {
    for (const float & tol : parameters_.tolerance.value()) {
      if (tol > 0) {
        toleranceVec_.push_back(tol);
      } else {
          throw eckit::UserError(
             "The tolerance value must be greater than zero.", Here());
      }
    }

    if (testVars_.nvars() != toleranceVec_.size()) {
      throw eckit::UserError(
        "The number of variables and the number of tolerances must be the same.", Here());
    }
  }

  timeToleranceSec_ = parameters_.analysisTimeTolerance.value().toSeconds();
  minSpacingSec_ = parameters_.minSpacing.value().toSeconds();
  const std::string selection = parameters_.equidistantTimeSelection.value();
  if (selection == "after") {
    preferAfter_ = true;
  } else if (selection == "before") {
    preferAfter_ = false;
  } else {
    throw eckit::UserError(
        "equidistantTimeSelection must be 'after' or 'before', got '" + selection + "'.",
          Here());
  }
  oops::Log::debug() << "DuplicateThinning: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateThinning::applyFilter(const std::vector<bool> & apply,
                                    const Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "DuplicateThinning applyFilter start" << std::endl;

  const util::DateTime missingDateTime = util::missingValue<util::DateTime>();
  const float missingFloat = util::missingValue<float>();
  const int missingInt = util::missingValue<int>();

  // Create ObsAccessor to get observations held in multiple MPI rank
  ObsAccessor obsAccessor = createObsAccessor();
  const size_t totalNumObs = obsAccessor.totalNumObservations();
  util::DateTime analysisTime = (parameters_.analysisTime.value() != boost::none)
                                ? *parameters_.analysisTime.value()
                                : obsdb_.windowStart();
  std::vector<bool> isThinned(totalNumObs, false);

  const std::vector<size_t> validObsIds
                  = obsAccessor.getValidObservationIds(apply);
  // Create RecursiveSplitter to partition array into groups of elements.
  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  std::vector<util::DateTime> dateTime = obsAccessor.getDateTimeVariableFromObsSpace(
                                          "MetaData", "dateTime");

  // Observations are grouped by the specified variable(s) and split.
  for (size_t i = 0; i < testVars_.nvars(); ++i) {
    if (obsdb_.dtype(testVars_[i].group(), testVars_[i].variable()) == ioda::ObsDtype::Float) {
      // 1. Read data values from ObsSpace
      std::vector<float> obsFloatCategories = obsAccessor.getFloatVariableFromObsSpace(
                      testVars_[i].group(), testVars_[i].variable());
      ASSERT(totalNumObs == obsFloatCategories.size());

      // 2. Configure tolerance
      auto binVal = [&](float v) -> int {
        if (v == missingFloat || std::isnan(v)) {
          return missingInt;
        }
        return static_cast<int>(std::lround(v / toleranceVec_[i]));
      };

      // 3. Project onto valid obs and bin the observations.
      std::vector<int> binnedCategories(validObsIds.size());
      for (size_t pos = 0; pos < validObsIds.size(); ++pos) {
        const size_t obsId = validObsIds[pos];
        float value = obsFloatCategories[obsId];
        binnedCategories[pos] = binVal(value);
      }
      // 4. Group observations by binned category and update splitter.
      splitter.groupBy(binnedCategories);
    } else if (obsdb_.dtype(testVars_[i].group(), testVars_[i].variable())
                == ioda::ObsDtype::Integer) {
      const std::vector<int> obsIntCategories = obsAccessor.getIntVariableFromObsSpace(
                      testVars_[i].group(), testVars_[i].variable());
      ASSERT(totalNumObs == obsIntCategories.size());
      splitter.groupBy(obsIntCategories);
    } else if (obsdb_.dtype(testVars_[i].group(), testVars_[i].variable())
                == ioda::ObsDtype::String) {
      const std::vector<std::string> obsStrCategories = obsAccessor.getStringVariableFromObsSpace(
                      testVars_[i].group(), testVars_[i].variable());
      ASSERT(totalNumObs == obsStrCategories.size());
      splitter.groupBy(obsStrCategories);
  } else {
      throw eckit::UserError("DuplicateThinning don't support dateTime variables"
           " for identifying duplicates.", Here());
    }
  }  // end for all variables

  // Sort groups by time and thin observation per time tolerance.
  splitter.sortGroupsBy([&dateTime, &validObsIds](size_t obsIndex)
                          { return dateTime[validObsIds[obsIndex]]; });

  for (RecursiveSplitter::Group categoryGroup : splitter.multiElementGroups()) {
    // First pass: apply the analysis-time-tolerance. skip and collect surviving obs IDs
    // for the second pass.
    std::vector<size_t> groupObsIds;
    groupObsIds.reserve(std::distance(categoryGroup.begin(), categoryGroup.end()));
    for (size_t pos : categoryGroup) {
      const size_t obsId = validObsIds[pos];
      // Obs with missing dateTime is not thinned.
      if (timeToleranceSec_ > 0 && dateTime[obsId] != missingDateTime) {
        const int64_t distFromSeed =
            std::abs((dateTime[obsId] - analysisTime).toSeconds());
        if (distFromSeed > timeToleranceSec_) {
          // Obs outside the analysis-time tolerance, thin it.
           isThinned[obsId] = true;
          continue;
        }
      }
      // Obs for second pass
      groupObsIds.push_back(obsId);
    }

    // Second pass: enforce minimum time spacing between retained obs inside group.
    // anything within `min_spacing` is thinned.
    // When two obs are equally close to analysis-time, prefered obs is chosen
    // using equidistant_time_selection.
    const bool preferAfter = preferAfter_;
    if (groupObsIds.size() > 1) {
      std::stable_sort(groupObsIds.begin(), groupObsIds.end(),
          [&dateTime, &analysisTime, &missingDateTime, preferAfter](size_t a, size_t b) {
            const int64_t da = (dateTime[a] == missingDateTime)
                ? std::numeric_limits<int64_t>::max()
                : std::abs((dateTime[a] - analysisTime).toSeconds());
            const int64_t db = (dateTime[b] == missingDateTime)
                ? std::numeric_limits<int64_t>::max()
                : std::abs((dateTime[b] - analysisTime).toSeconds());
            if (da != db) return da < db;  // choose closer value
            return preferAfter ? (dateTime[a] > dateTime[b])  // choose later value
                               : (dateTime[a] < dateTime[b]);  // choose earlier value
          });

      std::vector<util::DateTime> keptTime;
      keptTime.reserve(groupObsIds.size());
      for (size_t obsId : groupObsIds) {
        if (dateTime[obsId] == missingDateTime) {
          // Missing time — keep.
          continue;
        }
        bool tooClose = false;
        for (const util::DateTime & t : keptTime) {
          if (std::abs((dateTime[obsId] - t).toSeconds()) < minSpacingSec_) {
            tooClose = true;
            break;
          }
        }
        if (tooClose) {
          // Close to analysis time, thin the observation.
          isThinned[obsId] = true;
        } else {
          keptTime.push_back(dateTime[obsId]);
        }
      }
    }
  }

  // Flag any observations that have been thinned.
  obsAccessor.flagRejectedObservations(isThinned, flagged);
  oops::Log::trace() << "DuplicateThinning applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateThinning::print(std::ostream & os) const {
  os << "DuplicateThinning: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------
ObsAccessor DuplicateThinning::createObsAccessor() const {
  oops::Log::trace() << "DuplicateThinning createObsAccessor start" << std::endl;
  return ObsAccessor::toAllObservations(obsdb_);
}

}  // namespace ufo
