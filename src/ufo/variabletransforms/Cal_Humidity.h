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

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

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
* See VariableTransformsParameters for filter setup.
*/
class Cal_RelativeHumidity : public TransformBase {
 public:
  Cal_RelativeHumidity(const VariableTransformsParameters &options,
                       ioda::ObsSpace &os,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::vector<bool> &apply);
  // Run variable conversion
  void runTransform() override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT();
  void methodUKMO();
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
* See VariableTransformsParameters for filter setup.
*/
class Cal_SpecificHumidity : public TransformBase {
 public:
  Cal_SpecificHumidity(const VariableTransformsParameters &options,
                       ioda::ObsSpace &os,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::vector<bool> &apply);
  // Run check
  void runTransform() override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT();
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_HUMIDITY_H_
