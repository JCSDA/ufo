/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsDiagnostics.h"

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "ufo/SampledLocations.h"

#include "ioda/ObsSpace.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDiagnostics::ObsDiagnostics(const ioda::ObsSpace & os, const Locations_ & locs,
                               const oops::ObsVariables & vars)
  : obsdb_(os), gdiags_(locs, oops::Variables(vars.variables()))
{}

// -----------------------------------------------------------------------------

ObsDiagnostics::ObsDiagnostics(const eckit::Configuration & config, const ioda::ObsSpace & os,
                               const oops::ObsVariables & vars)
  : obsdb_(os), gdiags_(config, os, oops::Variables(vars.variables()))
{}

// -----------------------------------------------------------------------------

void ObsDiagnostics::allocate(const int nlev, const oops::ObsVariables & vars) {
  gdiags_.allocate(nlev, oops::Variables(vars.variables()));
}

// -----------------------------------------------------------------------------

void ObsDiagnostics::save(const std::vector<double> & vals,
                          const std::string & var,
                          const int lev) {
  gdiags_.putAtLevel(vals, oops::Variable(var), lev);
}

// -----------------------------------------------------------------------------

size_t ObsDiagnostics::nlevs(const std::string & var) const {
  return gdiags_.nlevs(oops::Variable(var));
}

// -----------------------------------------------------------------------------

void ObsDiagnostics::print(std::ostream & os) const {
  os << "ObsDiagnostics not printing yet.";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
