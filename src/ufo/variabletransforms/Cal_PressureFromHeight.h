/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_PRESSUREFROMHEIGHT_H_
#define UFO_VARIABLETRANSFORMS_CAL_PRESSUREFROMHEIGHT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/*!
* \brief Derive Pressure from height for vertical profile (e.g. sonde report)
*
* \details This is especially need for radiosonde using a 3 09 055 bufr
* template.
* To date there is only the UKMO formulation available.
*
* See VariableTransformsParameters for filter setup.
*/
class Cal_PressureFromHeightForProfile : public TransformBase {
 public:
  Cal_PressureFromHeightForProfile(const VariableTransformsParameters &options,
                                   ioda::ObsSpace &os,
                                   const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run variable conversion
  void runTransform() override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodUKMO();
};

/*!
* \brief Converts heights to pressures using the ICAO atmosphere.
*
* \details To date there is only the UKMO formulation available.
*
* See VariableTransformsParameters for filter setup.
*/
class Cal_PressureFromHeightForICAO : public TransformBase {
 public:
  Cal_PressureFromHeightForICAO(const VariableTransformsParameters &options,
                                ioda::ObsSpace &os,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run check
  void runTransform() override;

 private:
  // list of specific implementation(s) - This is controled by "method"
  void methodUKMO();
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_PRESSUREFROMHEIGHT_H_
