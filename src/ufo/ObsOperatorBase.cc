/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperatorBase.h"

#include <vector>

#include "eckit/config/Configuration.h"
#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "ufo/Locations.h"

namespace ufo {

// -----------------------------------------------------------------------------

std::unique_ptr<Locations> ObsOperatorBase::locations() const {
  std::vector<float> lons(odb_.nlocs());
  std::vector<float> lats(odb_.nlocs());
  std::vector<util::DateTime> times(odb_.nlocs());
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "dateTime", times);
  return std::unique_ptr<Locations>(new Locations(lons, lats, times, odb_.distribution()));
}

// -----------------------------------------------------------------------------

oops::Variables ObsOperatorBase::simulatedVars() const {
  return odb_.obsvariables();
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
                                             const ObsOperatorParametersBase & params) {
  oops::Log::trace() << "ObsOperatorBase::create starting" << std::endl;
  const std::string &id = params.name.value().value();
  typename std::map<std::string, ObsOperatorFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in ufo::ObsOperatorFactory." << std::endl;
    ABORT("Element does not exist in ufo::ObsOperatorFactory.");
  }
  ObsOperatorBase * ptr = jloc->second->make(odb, params);
  oops::Log::trace() << "ObsOperatorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsOperatorParametersBase>
ObsOperatorFactory::createParameters(const std::string &name) {
  typename std::map<std::string, ObsOperatorFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ufo::ObsOperatorFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
