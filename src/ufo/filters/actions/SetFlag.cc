/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/SetFlag.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/actions/RecordActionUtils.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

namespace {

bool isIgnored(IgnoredObservations ignoredObservations, int qcFlag) {
  switch (ignoredObservations) {
    case IgnoredObservations::NONE:
      return false;
    case IgnoredObservations::REJECTED:
      return QCflags::isRejected(qcFlag);
    case IgnoredObservations::DEFECTIVE:
      return QCflags::isDefective(qcFlag);
  }
  throw std::logic_error("Unhandled IgnoredObservations value");
}

}  // namespace

// -----------------------------------------------------------------------------

constexpr char IgnoredObservationsParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<IgnoredObservations>
  IgnoredObservationsParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------

static FilterActionMaker<SetFlag<true>> setFlagMaker_("set");
static FilterActionMaker<SetFlag<false>> unsetFlagMaker_("unset");

// -----------------------------------------------------------------------------

template <bool value>
SetFlag<value>::SetFlag(const SetFlagParameters &parameters)
  : parameters_(parameters) {
  oops::Log::trace() << "SetFlag constructor" << std::endl;
}

// -----------------------------------------------------------------------------

template <bool value>
bool SetFlag<value>::canSetOrUnsetAtLocation(int qcFlag) const {
  return !isIgnored(parameters_.ignore.value(), qcFlag);
}

// -----------------------------------------------------------------------------

template <bool value>
void SetFlag<value>::apply(const Variables &vars,
                           const std::vector<std::vector<bool>> &flagged,
                           ObsFilterData &data,
                           int,
                           ioda::ObsDataVector<int> &qcFlags,
                           ioda::ObsDataVector<float> &) const {
  if (parameters_.setFlagsToObservationReport && !parameters_.setObservationReportFlags) {
    throw eckit::UserError
      ("The option 'set variable flags to observation report' cannot be set to true if "
       "'set observation report flags' is set to false.", Here());
  }
  oops::Log::trace() << "SetFlag apply start" << std::endl;
  const std::string group = "DiagnosticFlags/" + parameters_.flag.value();

  const size_t nlocs = data.nlocs();
  std::vector<DiagnosticFlag> diagnosticFlags(nlocs);
  std::vector<DiagnosticFlag> diagnosticFlagsObsRep;

  if (parameters_.setObservationReportFlags) {
    if (!data.obsspace().has(group, "observationReport"))
      throw eckit::UserError("Variable '" + group + "/observationReport does not exist yet. "
                             "It needs to be set up with the 'Create Diagnostic Flags' filter "
                             "prior to using the 'set' or 'unset' action.");
    diagnosticFlagsObsRep.resize(nlocs);
    // Retrieve the current values of the diagnostic flag attached to the observation report.
    data.get(ufo::Variable(group + "/observationReport"), diagnosticFlagsObsRep);
  }

  // Loop over all filter variables
  for (size_t ifiltervar = 0, nvars = vars.nvars(); ifiltervar < nvars; ++ifiltervar) {
    const std::string variableName = vars.variable(ifiltervar).variable();
    if (!data.obsspace().has(group, variableName))
      throw eckit::UserError("Variable '" + group + '/' + variableName + "' does not exist yet. "
                             "It needs to be set up with the 'Create Diagnostic Flags' filter "
                             "prior to using the 'set' or 'unset' action.");
    // Retrieve the current values of the diagnostic flag attached to the current filter variable
    data.get(ufo::Variable(group + "/" + variableName), diagnosticFlags);
    // QC flags of the current filter variable
    const ioda::ObsDataRow<int> &filterVarQcFlags = qcFlags[variableName];
    for (size_t iobs = 0; iobs < nlocs; ++iobs) {
      // Set/unset the diagnostic flag if the filter has flagged this observation and
      // the action hasn't been told to skip it
      if (flagged[ifiltervar][iobs] && canSetOrUnsetAtLocation(filterVarQcFlags[iobs])) {
        diagnosticFlags[iobs] = value;
        if (parameters_.setObservationReportFlags && diagnosticFlagsObsRep[iobs] != value)
          diagnosticFlagsObsRep[iobs] = value;
      }
    }
    // Save the modified values of the diagnostic flag to the ObsSpace.
    // Do not do this if the flag will later be set to the value of the observation report
    // diagnostic flag.
    if (!parameters_.setFlagsToObservationReport)
      data.obsspace().put_db(group, variableName, diagnosticFlags,
                             vars.variable(ifiltervar).dimList());
  }

  if (parameters_.setObservationReportFlags) {
    data.obsspace().put_db(group, "observationReport", diagnosticFlagsObsRep);
    if (parameters_.setFlagsToObservationReport) {
      for (size_t ifiltervar = 0, nvars = vars.nvars(); ifiltervar < nvars; ++ifiltervar) {
        const std::string variableName = vars.variable(ifiltervar).variable();
        data.obsspace().put_db(group, variableName, diagnosticFlagsObsRep,
                               vars.variable(ifiltervar).dimList());
      }
    }
  }
  oops::Log::trace() << "SetFlag apply complete" << std::endl;
}

// -----------------------------------------------------------------------------

template <bool value>
void SetFlag<value>::apply_to_record(
    const Variables &vars, const std::vector<std::vector<bool>> &flagged,
    ObsFilterData &data, int /*filterQCflag*/,
    ioda::ObsDataVector<int> &qcFlags,
    ioda::ObsDataVector<float> &obserr) const {
  oops::Log::trace() << "SetFlag apply_to_record start" << std::endl;

  // Get record information: a vector for each record containing its observation
  // indices.
  const std::vector<std::vector<std::size_t>> recordLocs =
      actions::recordLocationsOrThrow(data, value ? "set" : "unset");

  // Pre-compute the index of each filter variable in the full QC flags array.
  std::vector<std::size_t> allVarIndexes;
  allVarIndexes.reserve(vars.nvars());
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    allVarIndexes.push_back(
        qcFlags.varnames().find(vars.variable(ifiltervar).variable()));
  }

  // Expand the per-location flagged mask to a whole-record mask.
  // Record expansion can only be triggered by locations that apply() would
  // actually update.
  const auto isEligibleForRecordExpansion = [&](std::size_t ifiltervar, std::size_t jloc) {
    return canSetOrUnsetAtLocation(qcFlags[allVarIndexes[ifiltervar]][jloc]);
  };
  const std::vector<std::vector<bool>> expandedFlagged =
      actions::expandFlaggedToWholeRecord(flagged, recordLocs, isEligibleForRecordExpansion);

  // Reuse per-observation logic in apply() once the record mask has been
  // expanded.
  apply(vars, expandedFlagged, data, 0, qcFlags, obserr);

  oops::Log::trace() << "SetFlag apply_to_record complete" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
