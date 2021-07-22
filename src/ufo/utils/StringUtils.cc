/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/StringUtils.h"

#include <sstream>
#include <string>

#include "eckit/exception/Exceptions.h"

namespace ufo {

// -----------------------------------------------------------------------------
// For now support both ioda v1 an ioda v2 syntax for specifying the combined
// variable and group name:
//     ioda v1: variable@group
//     ioda v2: group/variable
//
// We will eventually be obsoleting the ioda v1 syntax. For ioda v1 syntax, only
// allow for one "@" separator. For ioda v2 syntax, since it allows for nested groups,
// allow for multiple "/" separators, and the variable name is after the last "/".
//
void splitVarGroup(const std::string & vargrp, std::string & var, std::string & grp) {
  const size_t at = vargrp.find("@");
  const size_t slash = vargrp.find_last_of("/");

  if (at != std::string::npos) {
    // ioda v1 syntax
    var = vargrp.substr(0, at);
    grp = vargrp.substr(at + 1, std::string::npos);
    const size_t no_at = grp.find("@");
    ASSERT(no_at == std::string::npos);
  } else if (slash != std::string::npos) {
    // ioda v2 syntax
    grp = vargrp.substr(0, slash);
    var = vargrp.substr(slash + 1, std::string::npos);
  } else {
    // vargrp has no "@" nor "/" then assume that vargrp is a variable
    // name (with no group specified).
    var = vargrp;
    grp = "";
  }
}

// -----------------------------------------------------------------------------

void splitInstSat(const std::string & instsat, std::string & inst, std::string & sat) {
  const size_t at = instsat.find("_");
  inst = instsat.substr(0, at);
  sat = "";
  if (at != std::string::npos) {
    sat = instsat.substr(at + 1, std::string::npos);
    const size_t no_at = sat.find("_");
    ASSERT(no_at == std::string::npos);
  }
}

// -----------------------------------------------------------------------------

bool isFloat(const std::string & str) {
  std::istringstream iss(str);
  float factor;
  iss >> factor;
  return (iss.eof() && !iss.fail());
}

// -----------------------------------------------------------------------------

bool readFloat(const std::string & str, float & num) {
  std::istringstream iss(str);
  iss >> num;
  return (iss.eof() && !iss.fail());
}

}  // namespace ufo
