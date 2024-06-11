/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_VARIABLENAMEMAP_H_
#define UFO_UTILS_VARIABLENAMEMAP_H_

#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {

class VariableNameParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableNameParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<std::string> alias{"alias", this};
};

class VariableMapParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableMapParameters, Parameters)
 public:
  oops::OptionalParameter<std::vector<VariableNameParameters>> variableMap{"variable maps", this};
};

/// \brief Class which is used to allow different aliases for variable names. If any variables
/// require an alias these should be specified in a yaml file. The map of aliases is read from
/// the yaml file and stored in the Aliases_ member variable.
///
/// The map can be used to convert any variable name to an alias. If an alias file is not provided,
/// of if an alias file is provided but a variable alias is not provided, the name conversion
/// returns the original variable name.

class VariableNameMap {
 public:
    explicit VariableNameMap(const boost::optional<std::string> &);
    ~VariableNameMap();

  oops::Variable convertName(const std::string &) const;
  oops::Variables convertName(const oops::ObsVariables &) const;

 private:
  std::map<std::string, oops::Variable> Aliases_;
};

}  // namespace ufo

#endif  // UFO_UTILS_VARIABLENAMEMAP_H_
