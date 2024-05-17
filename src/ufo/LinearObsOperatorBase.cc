/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/LinearObsOperatorBase.h"

#include <memory>

#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

oops::ObsVariables LinearObsOperatorBase::simulatedVars() const {
  return odb_.assimvariables();
}

// -----------------------------------------------------------------------------

LinearObsOperatorFactory::LinearObsOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::LinearObsOperatorFactory."
                       << std::endl;
    ABORT("Element already registered in ufo::LinearObsOperatorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

LinearObsOperatorBase * LinearObsOperatorFactory::create(
    const ioda::ObsSpace & odb, const ObsOperatorParametersBase & params) {
  oops::Log::trace() << "LinearObsOperatorBase::create starting" << std::endl;

  const std::string &id = params.name.value().value();

  typename std::map<std::string, LinearObsOperatorFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in ufo::LinearObsOperatorFactory." << std::endl;
    ABORT("Element does not exist in ufo::LinearObsOperatorFactory.");
  }
  LinearObsOperatorBase * ptr = jloc->second->make(odb, params);
  oops::Log::trace() << "LinearObsOperatorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsOperatorParametersBase>
LinearObsOperatorFactory::createParameters(const std::string &name) {
  typename std::map<std::string, LinearObsOperatorFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ufo::LinearObsOperatorFactory");
  }
  return it->second->makeParameters();
}
// -----------------------------------------------------------------------------

}  // namespace ufo
