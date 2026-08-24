/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/DataHandlerParameters.h"

#include <algorithm>
#include <cstddef>
#include <string>

namespace ufo {

bool DataHandlerParameters::getOptional(const std::string &groupname) const {
  bool optional = false;
  if (std::find(groups_optional.value().begin(), groups_optional.value().end(), groupname)
      != groups_optional.value().end())
    optional = true;
  return optional;
}

size_t DataHandlerParameters::getEntriesPerProfile(const std::string &groupname) const {
  size_t entriesPerProfile = 0;
  // Variables with one entry per profile.
  if (std::find(groups_singlevalue.value().begin(),
                groups_singlevalue.value().end(), groupname)
      != groups_singlevalue.value().end()) {
    entriesPerProfile = 1;
  } else if (std::find(groups_modellevels.value().begin(),
                       groups_modellevels.value().end(), groupname)
      != groups_modellevels.value().end()) {
    entriesPerProfile = ModParameters.numModelLevels();
  } else if (std::find(groups_modelrholevels.value().begin(),
                       groups_modelrholevels.value().end(), groupname)
      != groups_modelrholevels.value().end()) {
    entriesPerProfile = ModParameters.numModelLevels_rho();
  }
  return entriesPerProfile;
}

}  // namespace ufo
