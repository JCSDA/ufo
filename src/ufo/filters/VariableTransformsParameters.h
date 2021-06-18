/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMSPARAMETERS_H_
#define UFO_FILTERS_VARIABLETRANSFORMSPARAMETERS_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
class Configuration;
}

namespace ufo {

/// \brief Options controlling the operation of the variablestansform filter.
class VariableTransformsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableTransformsParameters, Parameters)

 public:  // variables
  //=== Generic parameters ===//

  /// List of possible transformatioins:
  ///   - Pressure from height:
  ///        - \e "PressureFromHeightForProfile":\n
  ///          Retrieve pressure from observation height.
  ///          This conversion require a vertical profile. See
  ///          Cal_PressureFromHeightForProfile for details.
  ///        - \e "PressureFromHeightForICAO": \n
  ///          Retrieve pressure from observation height using ICAO standard.
  ///          See Cal_PressureFromHeightForICAO for details.
  ///        - \e "WindSpeedAndDirection": \n
  ///          Retrieve wind speed and direction from the eastward(u), and
  ///          northward (v) wind components.
  ///        - \e "WindComponents": \n
  ///          Retrieve the eastward(u), and northward (v) wind components from
  ///          wind speed and direction
  ///        - \e "SpecificHumidity": \n
  ///          Retrieve the specific humidity from relative humidity
  ///        - \e "RelativeHumidity": \n
  ///          Retrieve the relative humidity from specific humidity
  oops::RequiredParameter<std::vector<std::string>> Transform{"Transform",
                                                              this};

  /// Method used for calculation [Optional]:
  /// Related to Met Center - See ReadTheDoc for more details.
  oops::Parameter<std::string> Method{"Method", "default", this};

  /// Formulation possible [Optional]:
  /// By default \e formulation is set to \e the Method.
  /// See ReadTheDoc for more details
  oops::Parameter<std::string> Formulation{"Formulation", "", this};

  /// Should we use only the valid data? [Optional]:
  /// By default \e UseValidDataOnly is set to \e true.
  /// See ReadTheDoc for more details
  oops::Parameter<bool> UseValidDataOnly{"UseValidDataOnly", true, this};

  /// Should we allow super-saturated relative humidity? [Optional]:
  /// By default \e AllowSuperSaturation is set to \e false.
  /// See ReadTheDoc for more details
  oops::Parameter<bool> AllowSuperSaturation{"AllowSuperSaturation", false, this};
};
}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMSPARAMETERS_H_

