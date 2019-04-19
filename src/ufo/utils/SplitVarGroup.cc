/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/SplitVarGroup.h"

#include <string>

#include "eckit/exception/Exceptions.h"

namespace ufo {

// -----------------------------------------------------------------------------

void splitVarGroup(const std::string & vargrp, std::string & var, std::string & grp) {
  const size_t at = vargrp.find("@");
  var = vargrp.substr(0, at);
  if (at != std::string::npos) {
    grp = vargrp.substr(at + 1, std::string::npos);
    const size_t no_at = grp.find("@");
    ASSERT(no_at == std::string::npos);
  }
}

// -----------------------------------------------------------------------------
}  // namespace ufo
