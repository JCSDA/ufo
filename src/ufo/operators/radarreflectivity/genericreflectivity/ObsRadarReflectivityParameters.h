/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYPARAMETERS_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ReflectivityAlgorithmBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ReflectivityAlgorithmParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ReflectivityAlgorithmParametersWrapper, oops::Parameters)
 public:
  /// Name of the reflectivity algorithm.
  /// Valid names are specified using a `ReflectivityAlgorithmMaker` in subclasses of
  /// ReflectivityAlgorithmBase.
  oops::RequiredPolymorphicParameter<ReflectivityAlgorithmParametersBase,
    ReflectivityAlgorithmFactory>
    reflectivityAlgorithmName{"name", this};
};

/// Configuration options recognised by the RadarReflectivity operator.
class ObsRadarReflectivityParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadarReflectivityParameters, ObsOperatorParametersBase)

 public:
  /// Parameter that contains details of the reflectivity algorithm to use.
  oops::RequiredParameter<ReflectivityAlgorithmParametersWrapper>
    reflectivityAlgorithmParameters{"algorithm", this};

  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this operator is used as a
  /// component of the Composite operator. If `variables` is not set, the operator will simulate
  /// all ObsSpace variables. Please see the documentation of the Composite operator for further
  /// details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
     "variables",
     "List of variables to be simulated",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYPARAMETERS_H_
