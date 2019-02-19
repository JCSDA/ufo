/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperatorBase.h"

#include "eckit/config/Configuration.h"
#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "ufo/Locations.h"

namespace ufo {

// -----------------------------------------------------------------------------

Locations * ObsOperatorBase::locations(const util::DateTime & t1,
                                       const util::DateTime & t2) const {
  return new Locations(odb_, t1, t2);
}

// -----------------------------------------------------------------------------

ObsOperatorFactory::ObsOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::ObsOperatorFactory." << std::endl;
    ABORT("Element already registered in ufo::ObsOperatorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsOperatorBase * ObsOperatorFactory::create(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsOperatorBase::create starting" << std::endl;
  const std::string id = conf.getString("ObsType");
  typename std::map<std::string, ObsOperatorFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in ufo::ObsOperatorFactory." << std::endl;
    ABORT("Element does not exist in ufo::ObsOperatorFactory.");
  }
  ObsOperatorBase * ptr = jloc->second->make(odb, conf);
  oops::Log::trace() << "ObsOperatorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
