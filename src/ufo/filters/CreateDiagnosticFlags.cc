/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/CreateDiagnosticFlags.h"

#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/base/ObsVariables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/DiagnosticFlag.h"

namespace ufo {

// -----------------------------------------------------------------------------

CreateDiagnosticFlags::CreateDiagnosticFlags(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                                             std::shared_ptr<ioda::ObsDataVector<int>> qcflags,
                                             std::shared_ptr<ioda::ObsDataVector<float>> obserr)
  : ObsProcessorBase(obsdb, parameters.deferToPost, std::move(qcflags), std::move(obserr)),
    parameters_(parameters)
{
  oops::Log::trace() << "CreateDiagnosticFlags constructor starts" << std::endl;
  if (parameters.filterVariables.value().is_initialized()) {
    const oops::ObsVariables &observedVariables = obsdb.obsvariables();
    const oops::ObsVariables filterVariables = getFilterVariables();
    for (size_t i = 0; i < filterVariables.size(); ++i)
      if (!observedVariables.has(filterVariables[i])) {
        throw eckit::UserError("Filter variable '" + filterVariables[i] +
                               "' is not an observed variable", Here());
      }
  }
  oops::Log::trace() << "CreateDiagnosticFlags constructed" << std::endl;
}

// -----------------------------------------------------------------------------

CreateDiagnosticFlags::~CreateDiagnosticFlags() {
  oops::Log::trace() << "CreateDiagnosticFlags destructed" << std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables CreateDiagnosticFlags::getFilterVariables() const {
  if (parameters_.filterVariables.value().is_initialized()) {
    ufo::Variables vars;
    for (const Variable &var : *parameters_.filterVariables.value())
      vars += var;
    return vars.toOopsObsVariables();
  } else {
    return obsdb_.obsvariables();
  }
}

// -----------------------------------------------------------------------------

void CreateDiagnosticFlags::doFilter() {
  oops::Log::trace() << "CreateDiagnosticFlags doFilter starts" << std::endl;

  const oops::ObsVariables filterVars = getFilterVariables();
  // Loop over the names of flags we've been asked to create.
  for (const DiagnosticFlagParameters &params : parameters_.flags.value()) {
    const std::string &flagName = params.name;
    // Loop over observed variables.
    for (size_t jv = 0; jv < filterVars.size(); ++jv) {
      const std::string &varName = filterVars[jv];
      createFlag(flagName, varName, params.forceReinitialization, params.initialValue);
    }
    // Create observation report flags if required.
    if (parameters_.observationReportFlags) {
      const std::string varName = "observationReport";
      createFlag(flagName, varName, params.forceReinitialization, params.initialValue);
    }
  }

  oops::Log::trace() << "CreateDiagnosticFlags doFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void CreateDiagnosticFlags::print(std::ostream & os) const {
  os << "Create Diagnostic Flags: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

void CreateDiagnosticFlags::createFlag(const std::string & flagName, const std::string & varName,
                                       bool forceReinitialization, bool initialValue) const {
  bool createOrReinitializeVar = true;
  if (obsdb_.has("DiagnosticFlags/" + flagName, varName)) {
    if (!forceReinitialization)
      createOrReinitializeVar = false;
  }
  if (createOrReinitializeVar) {
    const std::string groupName = "DiagnosticFlags/" + flagName;
    // We need to create or reinitialize the Boolean variable corresponding to the
    // current flag and observed variable.
    obsdb_.put_db(groupName, varName,
                  std::vector<DiagnosticFlag>(obsdb_.nlocs(), initialValue));
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
