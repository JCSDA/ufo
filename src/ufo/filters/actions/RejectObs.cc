/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/RejectObs.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/actions/RecordActionUtils.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<RejectObs> makerRejectObs_("reject");

// -----------------------------------------------------------------------------

RejectObs::RejectObs(const RejectObsParameters &parameters)
  : allvars_(), parameters_(parameters) {
  oops::Log::trace() << "RejectObs constructor" << std::endl;
}

// -----------------------------------------------------------------------------

bool RejectObs::canRejectAtLocation(int currentFlag) {
  return currentFlag == QCflags::pass;
}

// -----------------------------------------------------------------------------

void RejectObs::apply(const Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      ObsFilterData &,
                      int filterQCflag,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  oops::Log::trace() << "RejectObs apply start" << std::endl;
  for (size_t jv = 0; jv < vars.nvars(); ++jv) {
    size_t iv = flags.varnames().find(vars.variable(jv).variable());
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[jv][jobs] && canRejectAtLocation(flags[iv][jobs]))
        flags[iv][jobs] = filterQCflag;
    }
  }
  oops::Log::trace() << "RejectObs apply complete" << std::endl;
}

// -----------------------------------------------------------------------------

void RejectObs::apply_to_record(const Variables &vars,
                                const std::vector<std::vector<bool>> &flagged,
                                ObsFilterData &data, int filterQCflag,
                                ioda::ObsDataVector<int> &flags,
                                ioda::ObsDataVector<float> &obserr) const {
  oops::Log::trace() << "RejectObs apply_to_record start" << std::endl;

  // Get record information: a vector for each record containing its observation indices.
  const std::vector<std::vector<std::size_t>> recordLocs =
      actions::recordLocationsOrThrow(data, "reject");

  // Pre-compute the index of each filter variable in the full QC flags array.
  // This avoids repeated lookups inside the predicate lambda.
  std::vector<std::size_t> allVarIndexes;
  allVarIndexes.reserve(vars.nvars());
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    allVarIndexes.push_back(flags.varnames().find(vars.variable(ifiltervar).variable()));
  }

  // Expand the per-location flagged mask to a whole-record mask.
  // Record expansion can only be triggered by locations that apply() would actually reject.
  const auto isEligibleForRecordExpansion = [&](std::size_t ifiltervar, std::size_t jloc) {
    return canRejectAtLocation(flags[allVarIndexes[ifiltervar]][jloc]);
  };
  const std::vector<std::vector<bool>> expandedFlagged =
      actions::expandFlaggedToWholeRecord(flagged, recordLocs, isEligibleForRecordExpansion);

  // Reuse per-observation logic in apply() once the record mask has been expanded.
  apply(vars, expandedFlagged, data, filterQCflag, flags, obserr);

  oops::Log::trace() << "RejectObs apply_to_record complete" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
