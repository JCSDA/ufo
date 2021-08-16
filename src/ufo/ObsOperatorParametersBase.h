/*
 * (C) Crown Copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSOPERATORPARAMETERSBASE_H_
#define UFO_OBSOPERATORPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {

// -----------------------------------------------------------------------------
/// \brief Base class of classes storing configuration parameters of specific observation operators
/// and linear observation operators.
class ObsOperatorParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsOperatorParametersBase, Parameters)

 public:
  /// \brief Observation operator type.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when observation operator parameters are deserialized into an
  /// ObsOperatorParametersWrapper and used by OperatorFactory to instantiate a operator whose
  /// type is determined at runtime), but not others (e.g. in tests written with a particular
  /// operator in mind). ObsOperatorParametersWrapper will throw an exception if this parameter
  /// is not provided.
  oops::OptionalParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ObsOperatorParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; operators using
/// GenericObsOperatorParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of ObsOperatorParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericObsOperatorParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericObsOperatorParameters, ObsOperatorParametersBase)
 public:
  oops::ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATORPARAMETERSBASE_H_
