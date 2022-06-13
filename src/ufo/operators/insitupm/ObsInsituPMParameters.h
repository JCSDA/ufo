/*
 * (C) Copyright 2021.
 *
 * This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_INSITUPM_OBSINSITUPMPARAMETERS_H_
#define UFO_OPERATORS_INSITUPM_OBSINSITUPMPARAMETERS_H_

#include <string>
#include <vector>

// List of Parameter classes
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the InsituPM operator.
class ObsInsituPMParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsInsituPMParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter <std::string> Model
    {"model",
     "Name of the model",
     this};

  oops::RequiredParameter <std::vector<std::string>> tracerGeovals
    {"tracer_geovals",
     "Names of model tracer variables",
     this};

  oops::Parameter <std::string> VertCoord
    {"vertical_coordinate",
     "Vertical coordinate for interpolation: "
     "height_asl or log_pressure",
     "height_asl", this};

  oops::OptionalParameter <std::vector<int>> tracerModesCMAQ
    {"tracer_modes_cmaq",
     "CMAQ PM modes: "
     "1-Aitken; 2-accumulation; 3-coarse",
     this};

  oops::Parameter <bool> UseScaleFacCMAQ
    {"use_scalefac_cmaq",
     "Use CMAQ scale factors or not: "
     "true/false - call PM25/PMtot routines",
     false, this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_INSITUPM_OBSINSITUPMPARAMETERS_H_
