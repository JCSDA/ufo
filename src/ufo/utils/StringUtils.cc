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

void splitVarGroup(const std::string & vargrp, std::string & var, std::string & grp) {
  const size_t at = vargrp.find("@");
  var = vargrp.substr(0, at);
  grp = "";
  if (at != std::string::npos) {
    grp = vargrp.substr(at + 1, std::string::npos);
    const size_t no_at = grp.find("@");
    ASSERT(no_at == std::string::npos);
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
