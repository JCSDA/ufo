/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "ufo/obsbias/ObsBiasBase.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasBase::ObsBiasBase(const eckit::Configuration & conf)
  : input_filename_(), output_filename_() {
  // Bias coefficients input file name
  if (conf.has("ObsBias.abias_in")) {
    input_filename_ = conf.getString("ObsBias.abias_in");
  }

  // Bias coefficients output file name
  if (conf.has("ObsBias.abias_out")) {
    output_filename_ = conf.getString("ObsBias.abias_out");
  }
}

// -----------------------------------------------------------------------------

ObsBiasFactory::ObsBiasFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::ObsBiasFactory." << std::endl;
    ABORT("Element already registered in ufo::ObsBiasFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsBiasBase * ObsBiasFactory::create(const eckit::Configuration & conf,
                                     const std::vector<std::string> & preds,
                                     const std::vector<int> & jobs) {
  oops::Log::trace() << "ObsBiasBase::create starting" << std::endl;
  if (conf.has("ObsBias")) {
    std::string id = "";
    id = conf.getString("ObsBias.name");
    typename std::map<std::string, ObsBiasFactory*>::iterator jloc = getMakers().find(id);
    if (jloc == getMakers().end()) {
      oops::Log::error() << id << " does not exist in ufo::ObsBiasFactory." << std::endl;
      ABORT("Element does not existed in ufo::ObsBiasFactory.");
    }
    ObsBiasBase * ptr = jloc->second->make(conf, preds, jobs);
    oops::Log::trace() << "ObsBiasBase::create done" << std::endl;
    return ptr;
  } else {
    return NULL;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
