/*
 * (C) Copyright 2020 Met Office UK
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ufo {

/// \brief Base class of classes storing parameters controlling specific observation filters.
class ObsFilterParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsFilterParametersBase, Parameters)
 public:
  /// \brief Observation filter type.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when observation filter parameters are deserialized into an
  /// ObsFilterParametersWrapper and used by FilterFactory to instantiate a filter whose type is
  /// determined at runtime), but not others (e.g. in tests written with a particular filter in
  /// mind). ObsFilterParametersWrapper will throw an exception if this parameter is not provided.
  oops::OptionalParameter<std::string> filter{"filter", this};
};

/// \brief A subclass of ObsFilterParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; filters using
/// GenericFilterParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of ObsFilterParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericObsFilterParameters : public ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericObsFilterParameters, ObsFilterParametersBase)
 public:
  oops::ConfigurationParameter config{this};
};

}  // namespace ufo
