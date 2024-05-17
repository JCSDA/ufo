/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/IodaGroupIndices.h"

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsGroup.h"
#include "oops/base/ObsVariables.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Returns indices of all elements in the range
/// [\p elements_to_look_for_begin, \p elements_to_look_for_end) in the \p all_elements vector,
/// in the same order that the former are in.
/// If \p throwexception is set, throws an exception if at least one of the elements looked for is
/// missing. Otherwise sets index as -1 for the missing elements.
template<typename T>
std::vector<int> getAllIndices(const std::vector<T> & all_elements,
                               typename std::vector<T>::const_iterator elements_to_look_for_begin,
                               typename std::vector<T>::const_iterator elements_to_look_for_end,
                               bool throwexception) {
  std::vector<int> result;
  for (typename std::vector<T>::const_iterator sought_it = elements_to_look_for_begin;
       sought_it != elements_to_look_for_end; ++sought_it) {
    const T &sought = *sought_it;
    const auto found_it = std::find(all_elements.begin(), all_elements.end(), sought);
    if (found_it != all_elements.end()) {
      result.push_back(std::distance(all_elements.begin(), found_it));
    } else {
      if (throwexception) {
        const std::string errormsg = "getAllIndices: Can't find element in the vector";
        throw eckit::BadParameter(errormsg, Here());
      } else {
        result.push_back(-1);
      }
    }
  }
  return result;
}

// -----------------------------------------------------------------------------
/// Returns indices of all elements in \p elements_to_look_for in the \p all_elements vector,
/// in the same order that \p elements_to_look_for are in.
/// If \p throwexception is set, throws an exception if at least one of the elements looked for is
/// missing. Otherwise sets index as -1 for the missing elements.
template<typename T>
std::vector<int> getAllIndices(const std::vector<T> & all_elements,
                               const std::vector<T> & elements_to_look_for,
                               bool throwexception) {
  return getAllIndices(all_elements, elements_to_look_for.begin(), elements_to_look_for.end(),
                       throwexception);
}

// -----------------------------------------------------------------------------
std::vector<int> getRequiredVariableIndices(const ioda::ObsGroup &obsgroup,
                 const std::string &varname,
                 typename std::vector<std::string>::const_iterator elements_to_look_for_begin,
                 typename std::vector<std::string>::const_iterator elements_to_look_for_end,
                 bool throwexception) {
  ioda::Variable var = obsgroup.vars.open(varname);
  std::vector<std::string> elements_in_group;
  var.read<std::string>(elements_in_group);
  return getAllIndices(elements_in_group, elements_to_look_for_begin, elements_to_look_for_end,
                       throwexception);
}

// -----------------------------------------------------------------------------

std::vector<int> getRequiredVarOrChannelIndices(const ioda::ObsGroup &obsgroup,
                                                const oops::ObsVariables &vars_to_look_for,
                                                bool throwexception) {
  if (vars_to_look_for.channels().empty()) {
    // Read all variables from the file into std vector
    // todo : pjn - need to change this to Variable at some point
    ioda::Variable variablesvar;
    if ( obsgroup.vars.exists("variables") ) {
       variablesvar = obsgroup.vars.open("variables");
    } else {
       if ( obsgroup.vars.exists("Variable") ) {
        variablesvar = obsgroup.vars.open("Variable");
       } else {
        const std::string msg("no variables or Variable found!");
        // todo pjn need to find the proper exception here
        throw eckit::Exception(msg, Here());
       }
    }
    std::vector<std::string> variables;
    variablesvar.read<std::string>(variables);

    // Find the indices of the ones we need
    return getAllIndices(variables, vars_to_look_for.variables(), throwexception);
  } else {
    // At present this function can search either for channel numbers or variable names,
    // but not both. So make sure there's only one multi-channel variable.
    ASSERT(vars_to_look_for.variables().size() == vars_to_look_for.channels().size());

    // Read all channels from the file into std vector
    // todo: pjn need to change this to Channel at some point
    ioda::Variable channelsvar;
    if ( obsgroup.vars.exists("channels") ) {
       channelsvar = obsgroup.vars.open("channels");
    } else {
       if ( obsgroup.vars.exists("Channel") ) {
        channelsvar = obsgroup.vars.open("Channel");
       } else {
        const std::string msg("no channels or Channel found!");
        // todo pjn need to find the proper exception here
        throw eckit::Exception(msg, Here());
       }
    }

    std::vector<int> channels;
    channelsvar.read<int>(channels);

    // Find the indices of the ones we need
    return getAllIndices(channels, vars_to_look_for.channels(), throwexception);
  }
}

// -----------------------------------------------------------------------------

std::vector<int> getAllStrIndices(const std::vector<std::string> & all_elements,
                   typename std::vector<std::string>::const_iterator elements_to_look_for_begin,
                   typename std::vector<std::string>::const_iterator elements_to_look_for_end,
                   bool throwexception) {
  return getAllIndices(all_elements, elements_to_look_for_begin, elements_to_look_for_end,
                       throwexception);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
