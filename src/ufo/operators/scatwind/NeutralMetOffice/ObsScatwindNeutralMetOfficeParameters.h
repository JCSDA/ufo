/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICEPARAMETERS_H_
#define UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICEPARAMETERS_H_

#include "oops/util/parameters/Parameters.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the neutral wind operator.
class ObsScatwindNeutralMetOfficeParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsScatwindNeutralMetOfficeParameters, ObsOperatorParametersBase)

 public:
  oops::Parameter<bool> surfaceTypeCheck
    {"surface_type_check",
     "If true, check that the surface type is sea before calculating H(x)."
     "Otherwise no check is performed.",
     true,
     this};
  oops::Parameter<int> surfaceTypeSea
    {"surface_type_sea",
     "The integer value used to denote sea in the surface type check.",
     0,
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICEPARAMETERS_H_
