/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
#define UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/stringFunctions.h"
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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const ufo::Variable &value) {
    eckit::LocalConfiguration subConfig;
    subConfig.set("name", value.variable() + "@" + value.group());
    const std::vector<int> &channels = value.channels();
    if (!channels.empty()) {
      const std::string channelsAsString = util::stringfunctions::join(
            ",", channels.begin(), channels.end(), [](int n) { return std::to_string(n); });
      subConfig.set("channels", channelsAsString);
    }
    if (!value.options().keys().empty())
      subConfig.set("options", value.options());
    config.set(name, subConfig);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    ObjectJsonSchema nestedSchema({{"name", {{"type", "\"string\""}}},
                                   {"options", {{"type", "\"object\""}}},
                                   {"channels", {{"type", "[\"string\", \"integer\"]"}}}});
    return ObjectJsonSchema({{name, nestedSchema.toPropertyJsonSchema()}});
  }
};

}  // namespace oops

#endif  // UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
