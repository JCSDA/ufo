/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsDiagnostics.h"

#include <ostream>
#include <string>

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

void ObsDiagnostics::save(const std::string &) const {}

// -----------------------------------------------------------------------------

void ObsDiagnostics::print(std::ostream & os) const {
  os << "ObsDiagnostics not printing yet.";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
