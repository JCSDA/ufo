/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLECONVERSIONPARAMETERS_H_
#define UFO_FILTERS_VARIABLECONVERSIONPARAMETERS_H_

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

/// \brief Options controlling the operation of the VariableConversion filter.
class VariableConversionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableConversionParameters, Parameters)

 public:  // variables
  //=== Generic parameters ===//

  /// List of possible transformatioins:
  ///   - Pressure from height:
  ///         - \e "PressureFromHeightForProfile":\n
  ///           Retrieve pressure from observation height.
  ///           This conversion require a vertical profile. See
  /// Cal_PressureFromHeightForProfile for details.
  ///         - \e "PressureFromHeightForICAO": \n
  ///           Retrieve pressure from observation height using ICAO standard.
  ///           See Cal_PressureFromHeightForICAO for details.
  oops::RequiredParameter<std::vector<std::string>> Calculate{"Calculate",
                                                              this};

  /// Method used for calculation [Optional]:
  ///   - UKMO
  ///   - NCAR
  ///   - NOAA
  ///   - DEFAULT
  ///
  /// By default \e method is set to \e "DEFAULT"
  oops::Parameter<std::string> Method{"Method", "default", this};
};
}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLECONVERSIONPARAMETERS_H_
