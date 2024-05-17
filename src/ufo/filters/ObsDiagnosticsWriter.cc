/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsDiagnosticsWriter.h"

#include "oops/util/Logger.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDiagnosticsWriter::ObsDiagnosticsWriter(
                       ioda::ObsSpace &, const Parameters_ & params,
                       std::shared_ptr<ioda::ObsDataVector<int> >,
                       std::shared_ptr<ioda::ObsDataVector<float> >)
  : params_(params), extradiagvars_()
{
  oops::Log::trace() << "ObsDiagnosticsWriter contructor" << std::endl;
  // Identify diagnostics variables
  Variables diagvars;
  if (params.filterVariables.value() != boost::none) {
  // read filter variables
    for (const Variable &var : *params.filterVariables.value())
      diagvars += var;
    extradiagvars_ += Variables(diagvars, "ObsDiag").toOopsObsVariables();
  }
}

// -----------------------------------------------------------------------------

void ObsDiagnosticsWriter::print(std::ostream & os) const {
  os << "ObsDiagnosticsWriter: " << params_.toConfiguration();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
