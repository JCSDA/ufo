/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/AcceptObs.h"

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "ufo/filters/actions/RecordActionUtils.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<AcceptObs> acceptObsMaker_("accept");

// -----------------------------------------------------------------------------

AcceptObs::AcceptObs(const Parameters_ &)
  : allvars_() {
  oops::Log::trace() << "AcceptObs constructor" << std::endl;
}

// -----------------------------------------------------------------------------

bool AcceptObs::canAcceptAtLocation(int currentFlag) {
  return currentFlag != QCflags::missing &&
         currentFlag != QCflags::preQC &&
         currentFlag != QCflags::Hfailed;
}

// -----------------------------------------------------------------------------

void AcceptObs::apply(const Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      ObsFilterData &,
                      int /*filterQCflag*/,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  oops::Log::trace() << "AcceptObs apply start" << std::endl;
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    const size_t iallvar = flags.varnames().find(vars.variable(ifiltervar).variable());
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[ifiltervar][jobs]) {
        int &currentFlag = flags[iallvar][jobs];
        if (canAcceptAtLocation(currentFlag))
          currentFlag = QCflags::pass;
      }
    }
  }
  oops::Log::trace() << "AcceptObs apply complete" << std::endl;
}

// -----------------------------------------------------------------------------

void AcceptObs::apply_to_record(const Variables &vars,
                                const std::vector<std::vector<bool>> &flagged,
                                ObsFilterData &data, int /*filterQCflag*/,
                                ioda::ObsDataVector<int> &flags,
                                ioda::ObsDataVector<float> &obserr) const {
  oops::Log::trace() << "AcceptObs apply_to_record start" << std::endl;

  // Get record information: a vector for each record containing its observation
  // indices.
  const std::vector<std::vector<std::size_t>> recordLocs =
      actions::recordLocationsOrThrow(data, "accept");

  // Pre-compute the index of each filter variable in the full QC flags array.
  std::vector<std::size_t> allVarIndexes;
  allVarIndexes.reserve(vars.nvars());
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    allVarIndexes.push_back(
        flags.varnames().find(vars.variable(ifiltervar).variable()));
  }

  // Expand the per-location flagged mask to a whole-record mask.
  // Record expansion can only be triggered by locations that apply() would
  // actually accept.
  const auto isEligibleForRecordExpansion = [&](std::size_t ifiltervar, std::size_t jloc) {
    return canAcceptAtLocation(flags[allVarIndexes[ifiltervar]][jloc]);
  };
  const std::vector<std::vector<bool>> expandedFlagged =
      actions::expandFlaggedToWholeRecord(flagged, recordLocs, isEligibleForRecordExpansion);

  // Reuse per-observation logic in apply() once the record mask has been
  // expanded. apply() handles the actual QC flag updates, so we avoid
  // duplicating that logic here.
  apply(vars, expandedFlagged, data, 0, flags, obserr);

  oops::Log::trace() << "AcceptObs apply_to_record complete" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
