/*
 * (C) Copyright 2020-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/AnalyticInitBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------

AnalyticInitFactory::AnalyticInitFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in analytic init factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<AnalyticInitBase> AnalyticInitFactory::create(const eckit::Configuration & config) {
  oops::Log::trace() << "AnalyticInitFactory::create starting" << std::endl;
  const std::string &id = config.getString("method");
  typename std::map<std::string, AnalyticInitFactory*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in analytic init factory.");
  }
  std::unique_ptr<AnalyticInitBase> ptr(jerr->second->make(config));
  oops::Log::trace() << "AnalyticInitFactory::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
