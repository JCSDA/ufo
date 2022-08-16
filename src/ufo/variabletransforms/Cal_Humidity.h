/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_HUMIDITY_H_
#define UFO_VARIABLETRANSFORMS_CAL_HUMIDITY_H_

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the variable transformation of specific humidity to
/// relative humidity
class Cal_HumidityParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_HumidityParameters, VariableTransformParametersBase);

 public:
  /// Should we allow super-saturated relative humidity? [Optional]:
  /// By default \e AllowSuperSaturation is set to \e false.
  /// See ReadTheDoc for more details
  oops::Parameter<bool> AllowSuperSaturation{"AllowSuperSaturation", false, this};
  oops::Parameter<std::string> SpecificHumidityVariable{"specific humidity variable",
                                                        "specificHumidity", this};
  oops::Parameter<std::string> PressureVariable{"pressure variable",
                                                "pressure", this};
  oops::Parameter<std::string> PressureAt2MVariable{"pressure at 2m variable",
                                                    "pressureAt2M", this};
  oops::Parameter<std::string> PressureGroupVariable{"pressure group variable",
                                                     "ObsValue", this};
  oops::Parameter<std::string> TemperatureVariable{"temperature variable",
                                                   "airTemperature", this};
  oops::Parameter<std::string> VirtualTempVariable{"virtual temperature variable",
                                                   "virtualTemperature", this};
  oops::Parameter<std::string> TemperatureAt2MVariable{"temperature at 2m variable",
                                                       "airTemperatureAt2M", this};
  oops::Parameter<std::string> RelativeHumidityVariable{"relative humidity variable",
                                                        "relativeHumidity", this};
  oops::Parameter<std::string> RelativeHumidityAt2MVariable{"relative humidity at 2m variable",
                                                            "relativeHumidityAt2M", this};
  oops::Parameter<std::string> WaterVaporMixingRatioVariable{"water vapor mixing ratio variable",
                                                             "waterVaporMixingRatio", this};
  oops::Parameter<std::string> DewPointTemperatureVariable{"dew point temperature variable",
                                                           "dewpointTemperature", this};
  oops::Parameter<std::string> DewPointTemperature2MVariable{"dew point temperature at 2m variable",
                                                             "dewpointTemperatureAt2M", this};
};

/*!
* \brief Relative Humidity filter
*
* Performs a variable conversion from specific_humidity, temperature, and
* pressure to relative humidity. The newly calculated variable is included in the same
* obs space.
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transform
*   Transform: ["RelativeHumidity"]
*   Formulation: Sonntag    # Using Sonntag formulation
* \endcode
*
* See Cal_RelativeHumidityParameters for filter setup.
*/
class Cal_RelativeHumidity : public TransformBase {
 public:
  typedef Cal_HumidityParameters Parameters_;

  Cal_RelativeHumidity(const Parameters_ &options,
                       const ObsFilterData &data,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  bool allowSuperSaturation_;
  std::string specifichumidityvariable_;
  std::string pressurevariable_;
  std::string pressureat2mvariable_;
  std::string pressuregroupvariable_;
  std::string temperaturevariable_;
  std::string temperatureat2mvariable_;
  std::string relativehumidityvariable_;
  std::string relativehumidityat2mvariable_;
  std::string watervapormixingratiovariable_;
  std::string dewpointtemperaturevariable_;
  std::string dewpointtemperatureat2mvariable_;
  std::string virtualtempvariable_;
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT(const std::vector<bool> &apply);
  void methodUKMOmixingratio(const std::vector<bool> &apply);
  void methodUKMO(const std::vector<bool> &apply);
};


/*!
* \brief Specific Humidity filter
*
* Performs a variable conversion from relative_humidity, temperature, and
* pressure to specific humidity. The newly calculated variable is included in the same
* obs space.
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transform
*   Transform: ["SpecificHumidity"]
*   Method: UKMO            # Using UKMO method and UKMO default formulation
*   Formulation: Sonntag    # Using Sonntag formulation
* \endcode
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_SpecificHumidity : public TransformBase {
 public:
  typedef Cal_HumidityParameters Parameters_;

  Cal_SpecificHumidity(const Parameters_ &options,
                       const ObsFilterData &data,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string specifichumidityvariable_;
  std::string pressurevariable_;
  std::string pressureat2mvariable_;
  std::string pressuregroupvariable_;
  std::string temperaturevariable_;
  std::string relativehumidityvariable_;
  std::string dewpointtemperaturevariable_;
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT(const std::vector<bool> &apply);
};

/*!
* \brief Virtual Temperature filter
*
* Performs a variable conversion from temperature and specific humidity to virtual temperature.
* The newly calculated variable is included in the same obs space.
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transform
*   Transform: "VirtualTemperature"
* \endcode
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_VirtualTemperature : public TransformBase {
 public:
  typedef Cal_HumidityParameters Parameters_;

  Cal_VirtualTemperature(const Parameters_ &options,
                         const ObsFilterData &data,
                         const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                         const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string specifichumidityvariable_;
  std::string temperaturevariable_;
  std::string virtualtempvariable_;
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT(const std::vector<bool> &apply);
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_HUMIDITY_H_
