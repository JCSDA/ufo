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

namespace ufo {

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
                           ioda::ObsDataVector<int> &,
                           ioda::ObsDataVector<float> &) const {
  const std::string group = "DiagnosticFlags/" + parameters_.flag.value();
  const size_t nlocs = data.nlocs();
  std::vector<DiagnosticFlag> flags(nlocs);
  for (size_t ivar = 0, nvars = vars.nvars(); ivar < nvars; ++ivar) {
    const std::string variableName = vars.variable(ivar).variable();
    if (!data.obsspace().has(group, variableName))
      throw eckit::UserError("Diagnostic flag '" + parameters_.flag.value() + "' does not exist "
                             "yet. It needs to be set up with the 'Create Diagnostic Flags' "
                             "filter prior to using the 'set' or 'unset' action.");
    data.get(ufo::Variable(group + "/" + variableName), flags);
    for (size_t iobs = 0; iobs < nlocs; ++iobs) {
      if (flagged[ivar][iobs]) {
        flags[iobs] = value;
      }
    }
    data.obsspace().put_db(group, variableName, flags);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
