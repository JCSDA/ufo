/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsFilterBase.h"


#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace ufo {

FilterFactory::FilterFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs filter factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsFilterBase>
FilterFactory::create(ioda::ObsSpace & os, const ObsFilterParametersBase & params,
                      ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr) {
  oops::Log::trace() << "FilterFactory::create starting" << std::endl;
  const std::string &id = params.filter.value().value();
  typename std::map<std::string, FilterFactory*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in obs filter factory." << std::endl;
    oops::Log::error() << "Obs Filter Factory has " << getMakers().size() <<
                          " elements:" << std::endl;
    for (typename std::map<std::string, FilterFactory*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       oops::Log::error() << "A " << jj->first << " Filter" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in obs filter factory.");
  }
  std::unique_ptr<ObsFilterBase> ptr(jloc->second->make(os, params, flags, obserr));
  oops::Log::trace() << "FilterFactory::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsFilterParametersBase>
FilterFactory::createParameters(const std::string &name) {
  typename std::map<std::string, FilterFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ObsFilterFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
