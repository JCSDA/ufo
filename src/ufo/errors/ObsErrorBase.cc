/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorBase.h"

#include <vector>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsErrorFactory::ObsErrorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::ObsErrorFactory." << std::endl;
    ABORT("Element already registered in ufo::ObsErrorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsErrorBase * ObsErrorFactory::create(const ObsErrorParametersBase &params,
                                       ioda::ObsSpace &odb) {
  oops::Log::trace() << "ObsErrorFactory::create starting" << std::endl;
  const std::string &id = params.model.value();
  typename std::map<std::string, ObsErrorFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in ufo::ObsErrorFactory." << std::endl;
    ABORT("Element does not exist in ufo::ObsErrorFactory.");
  }
  ObsErrorBase * ptr = jloc->second->make(params, odb);
  oops::Log::trace() << "ObsErrorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsErrorParametersBase>
ObsErrorFactory::createParameters(const std::string &name) {
  typename std::map<std::string, ObsErrorFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    oops::Log::error() << name << " does not exist in ufo::ObsErrorFactory." << std::endl;
    ABORT("Element does not exist in ufo::ObsErrorFactory.");
  }
  return it->second->makeParameters();
}
// -----------------------------------------------------------------------------

}  // namespace ufo
