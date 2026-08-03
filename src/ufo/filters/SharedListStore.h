/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_SHAREDLISTSTORE_H_
#define UFO_FILTERS_SHAREDLISTSTORE_H_

#include <map>
#include <string>

#include "eckit/config/LocalConfiguration.h"

namespace ufo {

/// \brief Singleton store for shared list configurations loaded from YAML files.
///
/// Lists are loaded lazily on first access and cached by file path.
/// Each MPI rank loads independently (read-only, no collision).
/// The store only caches the parsed YAML configuration — interpretation
/// (simple vs compound) is done by the filter at query time.
class SharedListStore {
 public:
  static SharedListStore & instance();

  /// Load a YAML file if not already loaded.
  void load(const std::string & filepath);

  /// Get the parsed configuration for a file.
  const eckit::LocalConfiguration & getConfig(const std::string & filepath) const;

  /// Check if a file has been loaded.
  bool isLoaded(const std::string & filepath) const;

 private:
  SharedListStore() = default;
  SharedListStore(const SharedListStore &) = delete;
  SharedListStore & operator=(const SharedListStore &) = delete;

  /// Cached configurations keyed by file path.
  std::map<std::string, eckit::LocalConfiguration> configs_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SHAREDLISTSTORE_H_
