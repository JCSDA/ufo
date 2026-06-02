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

#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the height to pressure conversion.
class Cal_PressureFromHeightParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_PressureFromHeightParameters, VariableTransformParametersBase);

 public:
  /// Height coordinate name.
  oops::RequiredParameter<std::string> HeightCoord{"height coordinate", this};

  /// Pressure coordinate name.
  oops::RequiredParameter<std::string> PressureCoord{"pressure coordinate", this};

  /// Pressure group name.
  oops::Parameter<std::string> PressureGroup{"pressure group", "ObsValue", this};
};


/*!
* \brief Derive Pressure from height for vertical profile (e.g. sonde report)
*
* \details This is especially need for radiosonde using a 3 09 055 bufr
* template.
* To date there is only the UKMO formulation available.
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_PressureFromHeightForProfile : public TransformBase {
 public:
  typedef Cal_PressureFromHeightParameters Parameters_;

  Cal_PressureFromHeightForProfile(const Parameters_ &options,
                                   const ObsFilterData &data,
                                   const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                   const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT(const std::vector<bool> &apply);

  /// Height coordinate name.
  std::string heightCoord_;

  /// Pressure coordinate name.
  std::string pressureCoord_;

  /// Pressure group name.
  std::string pressureGroup_;
};

/*!
* \brief Converts heights to pressures using the ICAO atmosphere.
*
* \details To date there is only the UKMO formulation available.
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_PressureFromHeightForICAO : public TransformBase {
 public:
  typedef Cal_PressureFromHeightParameters Parameters_;

  Cal_PressureFromHeightForICAO(const Parameters_ &options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;

 private:
  // list of specific implementation(s) - This is controlled by "method"
  void methodDEFAULT(const std::vector<bool> &apply);

  /// Height coordinate name.
  std::string heightCoord_;

  /// Pressure coordinate name.
  std::string pressureCoord_;

  /// Pressure group name.
  std::string pressureGroup_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_PRESSUREFROMHEIGHT_H_
