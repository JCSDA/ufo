/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_IODAGROUPINDICES_H_
#define UFO_UTILS_IODAGROUPINDICES_H_

#include <string>
#include <vector>

namespace ioda {
  class ObsGroup;
}

namespace oops {
  class Variables;
}

namespace ufo {

/// Return the vector of indices of all elements in the range
/// [\p elements_to_look_for_begin, \p elements_to_look_for_end) in the \p varname ioda::Variable
/// from \p obsgroup ioda::ObsGroup, in the same order that the former are in.
/// Throws an exception if at least one of elements looked for is missing.
std::vector<int> getRequiredVariableIndices(const ioda::ObsGroup &obsgroup,
                 const std::string &varname,
                 typename std::vector<std::string>::const_iterator elements_to_look_for_begin,
                 typename std::vector<std::string>::const_iterator elements_to_look_for_end);

/// Return the vector of indices of variables or channels \p vars_to_look_for in the
/// "variables" or "channels" ioda::Variable from \p obsgroup ioda::ObsGroup
std::vector<int> getRequiredVarOrChannelIndices(const ioda::ObsGroup &obsgroup,
                                                const oops::Variables &vars_to_look_for);

}  // namespace ufo

#endif  // UFO_UTILS_IODAGROUPINDICES_H_
