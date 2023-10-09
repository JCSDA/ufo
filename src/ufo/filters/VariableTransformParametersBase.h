/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMPARAMETERSBASE_H_
#define UFO_FILTERS_VARIABLETRANSFORMPARAMETERSBASE_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterParametersBase.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Abstract base class options controlling the operation of the variablestansform filter.
class VariableTransformParametersBase : public FilterParametersBase {
  OOPS_ABSTRACT_PARAMETERS(VariableTransformParametersBase, FilterParametersBase)

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
  oops::RequiredParameter<std::string> Transform{"Transform",
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

  /// Fill any missing entries of a vector in a Derived group (e.g. DerivedObsValue) with
  /// the non-missing entries of the vector in the equivalent original group (e.g. ObsValue).
  oops::Parameter<bool> FillMissingDerivedFromOriginal{"FillMissingDerivedFromOriginal",
                                                       false, this};

  /// Skip the variable transform calculation when there are no observations
  oops::Parameter<bool> SkipWhenNoObs{"SkipWhenNoObs", true, this};
};

/// \brief Concrete class containing the options specified by the VariableTransformParametersBase.
class GenericVariableTransformParameters : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericVariableTransformParameters, VariableTransformParametersBase)
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMPARAMETERSBASE_H_

