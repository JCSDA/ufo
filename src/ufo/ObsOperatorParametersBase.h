/*
 * (C) Crown Copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSOPERATORPARAMETERSBASE_H_
#define UFO_OBSOPERATORPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
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

  /// \brief Parameter specifying path to yaml file containing Observation to GeoVaL name mapping
  oops::OptionalParameter<std::string> AliasFile{"observation alias file", this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATORPARAMETERSBASE_H_
