/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/VariableNameMap.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/filesystem/PathName.h"
#include "oops/util/Logger.h"

namespace ufo {
VariableNameMap::VariableNameMap(const boost::optional<std::string> & aliasFile) {
  if (aliasFile.is_initialized()) {
    eckit::YAMLConfiguration conf(eckit::PathName(aliasFile.value()));
    VariableMapParameters mappingParams;
    mappingParams.validateAndDeserialize(conf);

    for (VariableNameParameters const& variableMap : *mappingParams.variableMap.value()) {
      Aliases_[variableMap.name] = variableMap.alias;
    }
  }
}

VariableNameMap::~VariableNameMap() {}

const std::string VariableNameMap::convertName(const std::string & name) {
  if (!Aliases_.empty()) {
    std::string alias;
    const auto & it = Aliases_.find(name);
    if (it == Aliases_.end()) {
      oops::Log::debug() << "Alias file supplied, but no alias found for variable " + name
                         << ". If your observation and GeoVaL names differ you will need to"
                            " include an alias for " << name << " in your name mapping file."
                         << std::endl;
      return name;
    } else {
      alias = it->second;
      oops::Log::debug() << "Variable " << name << " has been assigned alias "
                         << alias << std::endl;
    }
    return alias;
  } else {
    oops::Log::debug() << "No name mapping file supplied. If your observation and GeoVaL"
                          " names differ you will need to supply a name mapping file"
                          " containing name aliases." << std::endl;
    return name;
  }
}

oops::Variables VariableNameMap::convertName(const oops::Variables &vars) {
  std::vector<std::string> newvars;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    std::string alias = convertName(vars[jv]);
    newvars.push_back(alias);
  }
  oops::Variables aliasvars(newvars);
  return aliasvars;
}
}  // namespace ufo
