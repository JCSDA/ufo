/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/parameters/ParameterTraitsVariable.h"

#include <list>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/utils/StringTools.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/stringFunctions.h"
#include "ufo/filters/Variable.h"

namespace oops {

boost::optional<ufo::Variable> ParameterTraits<ufo::Variable>::get(
    util::CompositePath &path, const eckit::Configuration &config, const std::string& name) {
  if (config.has(name)) {
    std::list<std::string> messages;

    // Handle the following YAML structure:
    //
    // <name>:
    //   name: somevar@SomeGroup
    //   channels: ... # optional
    //   options:      # optional
    //     ...
    try {
      eckit::LocalConfiguration varConf(config, name);
      return ufo::Variable(varConf);
    } catch (eckit::Exception &e) {
      // The YAML doesn't have this structure.
      messages.emplace_back(e.what());
    }

    // Handle the following YAML structure:
    //
    // <name>: somevar@SomeGroup

    try {
      std::string varAndGroup = config.getString(name);
      return ufo::Variable(varAndGroup);
    } catch (eckit::Exception &e) {
      // The YAML doesn't have this structure.
      messages.emplace_back(e.what());
    }

    messages.push_front("The key '" + name +
                        "' is set to neither a string nor a map with the correct keys.");
    throw eckit::Exception(eckit::StringTools::join("\n", messages.begin(), messages.end()),
                           Here());
  } else {
    return boost::none;
  }
}

void ParameterTraits<ufo::Variable>::set(eckit::LocalConfiguration &config,
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

ObjectJsonSchema ParameterTraits<ufo::Variable>::jsonSchema(const std::string &name) {
  std::stringstream oneOf;
  {
    eckit::Channel ch;
    ch.setStream(oneOf);
    ch << "[\n";
    {
      eckit::AutoIndent indent(ch);
      ObjectJsonSchema simpleSchema = ParameterTraits<std::string>::jsonSchema("");
      ObjectJsonSchema completeSchema({{"name", {{"type", "\"string\""}}},
                                       {"options", {{"type", "\"object\""}}},
                                       {"channels", {{"type", "[\"string\", \"integer\"]"}}}});
      ch << toString(simpleSchema.properties().at("")) << ",\n";
      ch << completeSchema.toString() << '\n';
    }
    ch << "]";
  }

  return ObjectJsonSchema({{name, {{"oneOf", oneOf.str()}}}});
}

}  // namespace oops
