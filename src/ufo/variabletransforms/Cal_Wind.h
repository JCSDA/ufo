/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_WIND_H_
#define UFO_VARIABLETRANSFORMS_CAL_WIND_H_

#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the wind components transform.
class Cal_WindComponentsParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_WindComponentsParameters, VariableTransformParametersBase);

 public:
  /// Observation group name. Default is ObsValue.
  oops::Parameter<std::string> group{"group", "ObsValue", this};
  oops::Parameter<std::string> WindDirectionVariable{"wind direction variable",
                                                     "windDirection", this};
  oops::Parameter<std::string> WindSpeedVariable{"wind speed variable",
                                                "windSpeed", this};
};

/// Configuration parameters for the wind speed and direction transform.
class Cal_WindSpeedAndDirectionParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_WindSpeedAndDirectionParameters, VariableTransformParametersBase);

 public:
  /// Observation group name. Default is ObsValue.
  oops::Parameter<std::string> group{"group", "ObsValue", this};
  oops::Parameter<std::string> EastwardWindVariable{"eastward wind variable",
                                                    "windEastward", this};
  oops::Parameter<std::string> NorthwardWindVariable{"northward wind variable",
                                                     "windNorthward", this};
};

/*!
* \brief Wind Speed And Direction filter
*
* \details  Performs a variable conversion from the wind components, windEastward and
*  windNorthward, to windSpeed and windDirection. The newly calculated variables
*  are included in the same obs space. The filter does not have any configuration options.
///
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_WindSpeedAndDirection : public TransformBase {
 public:
  typedef Cal_WindSpeedAndDirectionParameters Parameters_;

  Cal_WindSpeedAndDirection(const Parameters_ &options,
                            const ObsFilterData &data,
                            const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                            const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string eastwardwindvariable_;
  std::string northwardwindvariable_;
  /// Group name.
  std::string group_;
};

/*!
* \brief Retrieve wind components.
*
* \details Performs a variable conversion from windSpeed and windDirection to
*  the wind components, windEastward and windNorthward. The newly calculated variables
*  are included in the same obs space. This filter supports the use of Channels as the
*  the second observed dimension. The variable conversion is performed for
*  each channel and so the output windEastward and windNorthward variables will have the
*  same dimensions as windSpeed and windDirection.
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_WindComponents : public TransformBase {
 public:
  typedef Cal_WindComponentsParameters Parameters_;

  Cal_WindComponents(const Parameters_ &options,
                     const ObsFilterData &data,
                     const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                     const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string windspeedvariable_;
  std::string winddirectionvariable_;
  /// Group name.
  std::string group_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_WIND_H_
