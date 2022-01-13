/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_POTENTIALTFROMT_H_
#define UFO_VARIABLETRANSFORMS_CAL_POTENTIALTFROMT_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

class Cal_PotentialTFromTParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_PotentialTFromTParameters, VariableTransformParametersBase);

 public:
  /// input pressure variable name
  oops::Parameter<std::string> PressureVariable{"pressure variable",
                                                "pressure", this};
  /// input pressure group name
  oops::Parameter<std::string> PressureGroup{"pressure group",
                                             "ObsValue", this};
  /// input temperature variable name
  oops::Parameter<std::string> TemperatureVariable{"temperature variable",
                                                   "airTemperature", this};
  /// potential temperature variable name
  oops::Parameter<std::string> PotentialTempVariable{"potential temperature variable",
                                                     "potentialTemperature", this};
};

/*!
* \brief Potential Temperature filter
*
* Performs a variable conversion from temperature, and
* pressure to potential temperature. The newly calculated variable is included in the same
* obs space.
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transforms
*   Transform: PotentialTFromT
* \endcode
*/

class Cal_PotentialTFromT : public TransformBase {
 public:
  typedef Cal_PotentialTFromTParameters Parameters_;

  Cal_PotentialTFromT(const Parameters_ &options,
                      const ObsFilterData &data,
                      const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                      const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  std::string pressurevariable_;
  std::string pressuregroup_;
  std::string temperaturevariable_;
  std::string potentialtempvariable_;
};
}  // namespace ufo
#endif  // UFO_VARIABLETRANSFORMS_CAL_POTENTIALTFROMT_H_
