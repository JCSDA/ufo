/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SharedListStore.h"

#include <string>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

SharedListStore & SharedListStore::instance() {
  static SharedListStore store;
  return store;
}

// -----------------------------------------------------------------------------

void SharedListStore::load(const std::string & filepath) {
  // Use basename as the cache key to avoid duplicates from different paths
  const std::string key = eckit::PathName(filepath).baseName();

  if (configs_.count(key) > 0) {
    oops::Log::debug() << "SharedListStore: using cached " << key << std::endl;
    return;
  }

  oops::Log::info() << "SharedListStore: loading " << filepath << std::endl;

  eckit::PathName path(filepath);
  eckit::YAMLConfiguration yamlConf(path);
  configs_[key] = eckit::LocalConfiguration(yamlConf);

  oops::Log::info() << "SharedListStore: loaded " << filepath
                    << " (key: " << key << ")" << std::endl;
}

// -----------------------------------------------------------------------------

const eckit::LocalConfiguration & SharedListStore::getConfig(
    const std::string & filepath) const {
  const std::string key = eckit::PathName(filepath).baseName();
  const auto it = configs_.find(key);
  if (it == configs_.end()) {
    throw eckit::UserError(
        "SharedListStore: file not loaded: " + filepath, Here());
  }
  return it->second;
}

// -----------------------------------------------------------------------------

bool SharedListStore::isLoaded(const std::string & filepath) const {
  const std::string key = eckit::PathName(filepath).baseName();
  return configs_.count(key) > 0;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
