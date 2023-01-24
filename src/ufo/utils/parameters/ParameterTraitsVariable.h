/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
#define UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_

#include <string>

#include "oops/util/parameters/ParameterTraits.h"
#include "ufo/filters/Variable.h"

/// \file ParameterTraitsVariable.h
/// This file needs to be included before any uses of Parameter<ufo::Variable> or
/// OptionalParameter<ufo::Variable>.

namespace oops {

/// Parameters representing variables can be set in YAML files in two ways:
///
/// * by setting a key to the variable name (a string), e.g.
///
///       selected variable: ObsValue/airTemperature
///
/// * by setting a key to a dictionary containing a key "name" set to the variable name
///   and optionally extra keys called "channels" and "options", e.g.
///
///       selected variable:
///         name: ObsValue/brightnessTemperature
///         channels: 1,5,7-10,13
///
///   or
///
///       selected variable:
///         name: ObsFunction/SCATRetMW
///         options:
///           scatret_ch238: 1
///           scatret_ch314: 2
///           scatret_ch890: 15
///           scatret_types: [ObsValue]
template <>
struct ParameterTraits<ufo::Variable> {
  static boost::optional<ufo::Variable> get(util::CompositePath &path,
                                            const eckit::Configuration &config,
                                            const std::string& name);

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const ufo::Variable &value);

  static ObjectJsonSchema jsonSchema(const std::string &name);
};

}  // namespace oops

#endif  // UFO_UTILS_PARAMETERS_PARAMETERTRAITSVARIABLE_H_
