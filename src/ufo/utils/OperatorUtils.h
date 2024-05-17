/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_OPERATORUTILS_H_
#define UFO_UTILS_OPERATORUTILS_H_

#include <vector>

#include <boost/optional.hpp>

#include "ufo/filters/Variable.h"

namespace eckit {
class Configuration;
}

namespace oops {
class ObsVariables;
}

namespace ioda {
class ObsSpace;
}

namespace ufo {
  /// Fill the list of variables to be simulated by an obs operator.
  ///
  /// \param listofvariables
  ///   An optional list of ufo::Variable objects that can be passed in via the configuration.
  ///   This controls which ObsSpace variables will be simulated. This option should
  ///   only be set if the operator is used as a component of the Composite operator.
  /// \param simulatedVariables
  ///   List of all simulated variables in the obs space associated with the obs operator.
  /// \param[out] operatorVariables
  ///   Set to the list of variables taken from \p listofvariables if this
  ///   option is present or to \p simulatedVariables if not.
  /// \param[out] operatorVariableIndices
  ///   Indices of the elements of \p simulatedVariables corresponding to the variables in
  ///   \p operatorVariables.
  void getOperatorVariables(const boost::optional<std::vector<ufo::Variable>> &listofvariables,
                            const oops::ObsVariables &simulatedVariables,
                            oops::ObsVariables &operatorVariables,
                            std::vector<int> &operatorVariableIndices);
}  // namespace ufo

#endif  // UFO_UTILS_OPERATORUTILS_H_
