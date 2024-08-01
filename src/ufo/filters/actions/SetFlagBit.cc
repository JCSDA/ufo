/*
 * (C) Copyright 2024, UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/SetFlagBit.h"

#include <string>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<SetFlagBit> setFlagMaker_("set flag bit");

// -----------------------------------------------------------------------------

SetFlagBit::SetFlagBit(const SetFlagBitParameters &parameters)
  : bitsetter_(1 << parameters.bit.value()) {
}

// -----------------------------------------------------------------------------

void SetFlagBit::apply(const Variables &vars,
                       const std::vector<std::vector<bool>> &flagged,
                       ObsFilterData &data,
                       int,
                       ioda::ObsDataVector<int> &qcFlags,
                       ioda::ObsDataVector<float> &) const {
  const std::string group = "DiagnosticFlags";
  const size_t nlocs = data.nlocs();
  const size_t nvars = vars.nvars();
  std::vector<int> diagnosticFlags(nlocs);
  // Loop over all filter variables
  for (size_t ifiltervar = 0; ifiltervar < nvars; ++ifiltervar) {
    const std::string variableName = vars.variable(ifiltervar).variable();
    if (!data.obsspace().has(group, variableName))
      throw eckit::UserError("Variable '" + group + '/' + variableName + "' does not exist yet. "
                             "It needs to be set up with the 'Create Diagnostic Flags' filter "
                             "prior to using the 'set flag bit' action.");
    // Retrieve the current values of the diagnostic flag attached to the current filter variable
    data.get(ufo::Variable(group + "/" + variableName), diagnosticFlags);
    for (size_t iobs = 0; iobs < nlocs; ++iobs) {
      // Set/unset the diagnostic flag bit if the filter has flagged this observation
      if (flagged[ifiltervar][iobs]) {
        diagnosticFlags[iobs] |= bitsetter_;
      }
    }
    data.obsspace().put_db(group, variableName, diagnosticFlags);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
