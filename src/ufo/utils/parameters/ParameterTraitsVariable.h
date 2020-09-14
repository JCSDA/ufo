/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
#define UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_

#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "ufo/filters/Variable.h"

/// \file ParameterTraitsVariable.h
/// This file needs to be included before any uses of Parameter<ufo::Variable> or
/// OptionalParameter<ufo::Variable>.

namespace oops {

template <>
struct ParameterTraits<ufo::Variable> {
  static boost::optional<ufo::Variable> get(util::CompositePath &path,
                                            const eckit::Configuration &config,
                                            const std::string& name) {
    if (config.has(name)) {
      eckit::LocalConfiguration varConf(config, name);
      if (!varConf.has("name")) {
        // TODO(wsmigaj): shouldn't ufo::Variable itself throw an exception if
        // the 'name' property is not specified?
        throw eckit::BadParameter(path.path() + ": No variable name specified", Here());
      }
      return ufo::Variable(varConf);
    } else {
      return boost::none;
    }
  }
};

}  // namespace oops

#endif  // UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
