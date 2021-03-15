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

PredictorBase::PredictorBase(const eckit::Configuration & conf, const oops::Variables & vars)
  : func_name_(conf.getString("predictor.name")),
    geovars_(), hdiags_(), vars_(vars) {
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

PredictorBase * PredictorFactory::create(const eckit::Configuration & conf,
                                         const oops::Variables & vars) {
  oops::Log::trace() << "PredictorBase::create starting" << std::endl;
  const std::string name = conf.getString("predictor.name");
  if (!predictorExists(name)) {
    oops::Log::error() << name << " does not exist in ufo::PredictorFactory."
                       << std::endl;
    ABORT("Element does not exist in ufo::PredictorFactory.");
  }
  typename std::map<std::string, PredictorFactory*>::iterator jloc =
           getMakers().find(name);
  PredictorBase * ptr = jloc->second->make(conf, vars);
  oops::Log::trace() << "PredictorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

bool PredictorFactory::predictorExists(const std::string & name) {
  return (getMakers().find(name) != getMakers().end());
}

// -----------------------------------------------------------------------------

}  // namespace ufo
