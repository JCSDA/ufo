/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <string>

#include "ioda/ObsSpace.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/OperatorUtils.h"

namespace ufo {

void getOperatorVariables(const boost::optional<std::vector<ufo::Variable>> &listofvariables,
                          const oops::ObsVariables &simulatedVariables,
                          oops::ObsVariables &operatorVariables,
                          std::vector<int> &operatorVariableIndices) {
  if (listofvariables.is_initialized() && (listofvariables.value().size() > 0)) {
    Variables optionalvariables(listofvariables.get());
    operatorVariables = optionalvariables.toOopsObsVariables();

    std::map<std::string, int> simulatedVariableIndices;
    for (size_t i = 0; i < simulatedVariables.size(); ++i)
      simulatedVariableIndices[simulatedVariables[i]] = i;

    for (size_t i = 0; i < operatorVariables.size(); ++i) {
      auto it = simulatedVariableIndices.find(operatorVariables[i]);
      if (it == simulatedVariableIndices.end())
        throw eckit::BadValue("Operator variable '" + operatorVariables[i] +
                              "' isn't one of the simulated variables", Here());
      operatorVariableIndices.push_back(it->second);
    }
  } else {
    operatorVariables = simulatedVariables;
    operatorVariableIndices.resize(operatorVariables.size());
    std::iota(operatorVariableIndices.begin(), operatorVariableIndices.end(), 0);
  }
}

}  // namespace ufo
