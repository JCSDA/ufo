/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/FilterActionBase.h"

#include <map>
#include <string>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterActionFactory::FilterActionFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::FilterActionFactory." << std::endl;
    ABORT("Element already registered in ufo::FilterActionFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

FilterActionBase * FilterActionFactory::create(const eckit::Configuration & conf) {
  oops::Log::trace() << "FilterActionBase::create starting" << std::endl;
  const std::string name = conf.getString("name", "reject");
  typename std::map<std::string, FilterActionFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    oops::Log::error() << name << " does not exist in ufo::FilterActionFactory." << std::endl;
    ABORT("Element does not exist in ufo::FilterActionFactory.");
  }
  FilterActionBase * ptr = jloc->second->make(conf);
  oops::Log::trace() << "FilterActionBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
