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

#include "oops/base/Variables.h"
#include "ufo/Locations.h"

#include "ioda/ObsSpace.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDiagnostics::ObsDiagnostics(const ioda::ObsSpace & os, const Locations & locs,
                               const oops::Variables & vars)
  : obsdb_(os), gdiags_(locs, vars)
{}

// -----------------------------------------------------------------------------

ObsDiagnostics::ObsDiagnostics(const eckit::Configuration & conf, const ioda::ObsSpace & os,
                               const oops::Variables & vars)
  : obsdb_(os), gdiags_(conf, os, vars)
{}

// -----------------------------------------------------------------------------

void ObsDiagnostics::save(const std::vector<double> & vals,
                          const std::string & var,
                          const int lev) {
  gdiags_.put(vals, var, lev);
}

// -----------------------------------------------------------------------------

size_t ObsDiagnostics::nlevs(const std::string & var) const {
  return gdiags_.nlevs(var);
}

// -----------------------------------------------------------------------------

void ObsDiagnostics::get(std::vector<float> & vals, const std::string & var) const {
  gdiags_.get(vals, var);
}

// -----------------------------------------------------------------------------

void ObsDiagnostics::get(std::vector<float> & vals, const std::string & var,
                         const int lev) const {
  gdiags_.get(vals, var, lev);
}

// -----------------------------------------------------------------------------

void ObsDiagnostics::print(std::ostream & os) const {
  os << "ObsDiagnostics not printing yet.";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
