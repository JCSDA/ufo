/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_OPERATORUTILS_H_
#define UFO_UTILS_OPERATORUTILS_H_

#include <vector>

namespace eckit {
class Configuration;
}

namespace oops {
class Variables;
}

namespace ioda {
class ObsSpace;
}

namespace ufo {
  /// Fill the list of variables to be simulated by an obs operator.
  ///
  /// \param conf
  ///   Configuration of the obs operator.
  /// \param simulatedVariables
  ///   List of all simulated variables in the obs space associated with the obs operator.
  /// \param[out] operatorVariables
  ///   Set to the list of variables taken from the `variables` option in \p conf if this
  ///   option is present or to \p simulatedVariables if not.
  /// \param[out] operatorVariableIndices
  ///   Indices of the elements of \p simulatedVariables corresponding to the variables in
  ///   \p operatorVariables.
  void getOperatorVariables(const eckit::Configuration &conf,
                            const oops::Variables &simulatedVariables,
                            oops::Variables &operatorVariables,
                            std::vector<int> &operatorVariableIndices);
}  // namespace ufo

#endif  // UFO_UTILS_OPERATORUTILS_H_
