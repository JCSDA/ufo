/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/PredictorBase.h"

#include <map>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

PredictorBase::PredictorBase(const PredictorParametersBase & parameters,
                             const oops::ObsVariables & vars)
  : vars_(vars), geovars_(), hdiags_(), func_name_(parameters.name) {
}

// -----------------------------------------------------------------------------

PredictorFactory::PredictorFactory(const std::string & name) {
  if (predictorExists(name)) {
    oops::Log::error() << name << " already registered in ufo::PredictorFactory."
                               << std::endl;
    ABORT("Element already registered in ufo::PredictorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<PredictorBase> PredictorFactory::create(const PredictorParametersBase & parameters,
                                                        const oops::ObsVariables & vars) {
  oops::Log::trace() << "PredictorBase::create starting" << std::endl;
  const std::string name = parameters.name;
  if (!predictorExists(name)) {
    oops::Log::error() << name << " does not exist in ufo::PredictorFactory."
                       << std::endl;
    ABORT("Element does not exist in ufo::PredictorFactory.");
  }
  typename std::map<std::string, PredictorFactory*>::iterator jloc =
           getMakers().find(name);
  std::unique_ptr<PredictorBase> ptr = jloc->second->make(parameters, vars);
  oops::Log::trace() << "PredictorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<PredictorParametersBase> PredictorFactory::createParameters(
    const std::string &name) {
  typename std::map<std::string, PredictorFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ufo::PredictorFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

bool PredictorFactory::predictorExists(const std::string & name) {
  return (getMakers().find(name) != getMakers().end());
}

// -----------------------------------------------------------------------------

}  // namespace ufo
