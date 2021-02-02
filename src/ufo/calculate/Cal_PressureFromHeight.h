/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CALCULATE_CAL_PRESSUREFROMHEIGHT_H_
#define UFO_CALCULATE_CAL_PRESSUREFROMHEIGHT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/calculate/CalculateBase.h"

namespace ufo {

/*!
* \brief Derive Pressure from height for vertical profile (e.g. sonde report)
*
* \details This is especially need for radiosonde using a 3 09 055 bufr
* template.
* To date there is only the UKMO formulation available.
*
* See VariableConversionParameters for filter setup.
*/
class Cal_PressureFromHeightForProfile : public CalculateBase {
 public:
  Cal_PressureFromHeightForProfile(const VariableConversionParameters &options,
                                   ioda::ObsSpace &os,
                                   const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run variable conversion
  void runCalculate() override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodUKMO();
};

/*!
* \brief Converts heights to pressures using the ICAO atmosphere.
*
* \details To date there is only the UKMO formulation available.
*
* See VariableConversionParameters for filter setup.
*/
class Cal_PressureFromHeightForICAO : public CalculateBase {
 public:
  Cal_PressureFromHeightForICAO(const VariableConversionParameters &options,
                                ioda::ObsSpace &os,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run check
  void runCalculate() override;

 private:
  // list of specific implementation(s) - This is controled by "method"
  void methodUKMO();
};
}  // namespace ufo

#endif  // UFO_CALCULATE_CAL_PRESSUREFROMHEIGHT_H_
