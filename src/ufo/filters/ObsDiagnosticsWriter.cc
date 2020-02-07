/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsDiagnosticsWriter.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDiagnosticsWriter::ObsDiagnosticsWriter(
                       ioda::ObsSpace &, const eckit::Configuration & config,
                       boost::shared_ptr<ioda::ObsDataVector<int> >,
                       boost::shared_ptr<ioda::ObsDataVector<float> >)
  : config_(config), extradiagvars_()
{
  oops::Log::trace() << "ObsDiagnosticsWriter contructor" << std::endl;
  if (config_.has("filter variables")) {
    Variables diagvars(config_.getSubConfigurations("filter variables"));
    extradiagvars_ += Variables(diagvars, "ObsDiag").toOopsVariables();
  }
}

// -----------------------------------------------------------------------------

void ObsDiagnosticsWriter::print(std::ostream & os) const {
  os << "ObsDiagnosticsWriter: " << config_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
