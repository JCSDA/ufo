/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/SetFlag.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

namespace {

bool returnFalse(int /*qcflag*/) {
  return false;
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
}

// -----------------------------------------------------------------------------

template <bool value>
void SetFlag<value>::apply(const Variables &vars,
                           const std::vector<std::vector<bool>> &flagged,
                           const ObsFilterData &data,
                           int,
                           ioda::ObsDataVector<int> &qcFlags,
                           ioda::ObsDataVector<float> &) const {
  const std::string group = "DiagnosticFlags/" + parameters_.flag.value();

  typedef bool (*Predicate)(int);
  // Pointer to a function taking a QC flag and returning true if observations with this QC flag
  // should be ignored and false otherwise
  Predicate isIgnored;
  switch (parameters_.ignore.value()) {
  case IgnoredObservations::NONE:
    isIgnored = &returnFalse;
    break;

  case IgnoredObservations::REJECTED:
    isIgnored = &QCflags::isRejected;
    break;

  case IgnoredObservations::DEFECTIVE:
    isIgnored = &QCflags::isDefective;
    break;
  }

  const size_t nlocs = data.nlocs();
  std::vector<DiagnosticFlag> diagnosticFlags(nlocs);

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
      if (flagged[ifiltervar][iobs] && !isIgnored(filterVarQcFlags[iobs])) {
        diagnosticFlags[iobs] = value;
      }
    }
    // Save the modified values of the diagnostic flag to the ObsSpace
    data.obsspace().put_db(group, variableName, diagnosticFlags);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
