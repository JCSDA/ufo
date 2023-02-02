/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_SATZENITHANGLECORRECTION_H_
#define UFO_VARIABLETRANSFORMS_CAL_SATZENITHANGLECORRECTION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the satellite zenith angle correction variable transformation
class Cal_SatZenithAngleCorrectionParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_SatZenithAngleCorrectionParameters, VariableTransformParametersBase);

 public:
  /// The variable that will be transformed.  This will take the same format as other variables e.g
  /// transform variable:
  ///   name: surface_emissivity@VarMetaData
  ///   channels: 1-20
  oops::RequiredParameter<Variable> transformVariable{"transform variable", this};

  /// The transform that is performed by this filter is of the form:
  /// var = var + coeffA * sza^exponentA + coeffB * sza^exponentB +
  ///       coeffC * sza^exponentC
  /// where sza is the sensor zenith angle and ^ represents raising to the power.
  /// The below variables are the necessary coefficients and exponents for this
  /// conversion.  If the transform variable has channels then the coeffs will
  /// need to be of the same size as the input channels.
  oops::Parameter<std::vector<float>> coeffA{"coefficient a", {1}, this};
  oops::Parameter<std::vector<float>> coeffB{"coefficient b", {1}, this};
  oops::Parameter<std::vector<float>> coeffC{"coefficient c", {1}, this};
  oops::Parameter<int> exponentA{"exponent a", 1, this};
  oops::Parameter<int> exponentB{"exponent b", 1, this};
  oops::Parameter<int> exponentC{"exponent c", 1, this};

  /// The minimum value that the output var can be. If the input var or sensor zenith angle
  /// are missing the value of var will be missing.
  oops::OptionalParameter<float> minvalue{"minimum value", this};

  /// The maximum value that the output var can be. If the input var or sensor zenith angle
  /// are missing the value of var will be missing.
  oops::OptionalParameter<float> maxvalue{"maximum value", this};
};

/*!
* \brief Sensor zenith angle empirical correction
*
* Performs a variable transform using a conversion of the type:
* var = var + coeffA * sza^exponentA + coeffB * sza^exponentB + coeffC * sza^exponentC
* where var is the variable to be transformed, sza is the sensor zenith angle.
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variable Transforms
*   Transform: SatZenithAngleCorrection
*   transform variable:
*     name: DerivedObsValue/emissivity
*     channels: 1,2
*   ## Test coeffs
*   coefficient a: [-3.60E-03, -2.38E-03]
*   coefficient b: [ 2.21E-05,  2.15E-05]
*   coefficient c: [-7.83E-09, -5.00E-09]
*   exponent a: 0
*   exponent b: 2
*   exponent c: 4
*   minimum value: 0.0
*   maximum value: 1.0
* \endcode
*
* See Cal_SatZenithAngleCorrectionParameters for filter setup.
*/

class Cal_SatZenithAngleCorrection : public TransformBase {
 public:
  typedef Cal_SatZenithAngleCorrectionParameters Parameters_;

  Cal_SatZenithAngleCorrection(const Parameters_ &options,
                               const ObsFilterData &data,
                               const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                               const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
  Variables requiredVariables() const override { return variables_; }
 private:
  Parameters_ parameters_;
  Variables variables_;
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_SATZENITHANGLECORRECTION_H_
