/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMSCALEDRADIANCE_H_
#define UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMSCALEDRADIANCE_H_

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

/// Configuration parameters for the scaled radiance to radiance variable transformation
class Cal_SatRadianceFromScaledRadianceParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_SatRadianceFromScaledRadianceParameters,
                           VariableTransformParametersBase);

 public:
  /// The variable that will be transformed.  This will take the same format as other variables e.g
  /// transform variable:
  ///   name: ObsValue/radiance
  ///   channels: 1-20
  oops::RequiredParameter<Variable> transformVariable{"transform from", this};

  /// Radiance scaling is done using scale factors and a range of channels with which this is
  /// applied. The following variables are required to specify how many scale factors and where
  /// they are in the obs space.  The names are the base name i.e. without numbers attached.
  /// e.g. MetData/scaleFactor and numScaleFactors = 2 will read in MetaData/scaleFactor1 and
  /// MetaData/scaleFactor2
  oops::RequiredParameter<size_t> numScaleFactors{"number of scale factors", this};
  oops::RequiredParameter<Variable> scalingVariable{"scale factor variable", this};
  oops::RequiredParameter<Variable> scalingStartChannel{"scale factor start", this};
  oops::RequiredParameter<Variable> scalingEndChannel{"scale factor end", this};

  /// By default (false) this option will read the radiance scaling factors and start and end
  /// indexes from a single array e.g. scalingVariable = channelScaleFactor
  /// channelScaleFactor - which is of size nfactors, where nfactors must equal numScaleFactors.
  /// If set to true the radiance scale factors are expected to be in arrays e.g.
  /// scalingVariable = channelScaleFactor and scalingStartChannel = channelScaleFactor then
  /// channelScaleFactor1, channelScaleFactor2, ..., channelScaleFactor(numScaleFactors) and
  /// startChannelScale1, startChannelScale2, ..., startChannelScale(numScaleFactors)
  /// where channelScaleFactor1 and startChannelScale1, for example, are of size nlocs.
  oops::Parameter<bool> getFactorsFromMultipleArrays{"get scaling factors from multiple arrays",
                                                     false, this};
};

/*!
* \brief Convert scaled radiance to radiance.
*
* Performs a conversion of a scaled radiance to an unscaled radiance.  This takes the form:
* radiance = (scaled radiance)x10^(-channelScaleFactor)
*
* Example using all the parameter options:
*
* \code{.yaml}
* obs filters:
* - filter: Variable Transforms
*   Transform: SatRadianceFromScaledRadiance
*   transform from:
*     name: ObsValue/scaledRadiance
*     channels: *all_channels
*   number of scale factors: 10
*   scale factor variable: MetaData/channelScaleFactor
*   scale factor start: MetaData/startChannel
*   scale factor end: MetaData/endChannel
*   get scaling factors from multiple arrays: true
* \endcode
*
* See Cal_SatRadianceFromScaledRadianceParameters for filter setup.
*/

class Cal_SatRadianceFromScaledRadiance : public TransformBase {
 public:
  typedef Cal_SatRadianceFromScaledRadianceParameters Parameters_;

  Cal_SatRadianceFromScaledRadiance(const Parameters_ &options,
                                 const ObsFilterData &data,
                                 const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                 const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
  Variables requiredVariables() const override { return variables_; }
 private:
  Parameters_ parameters_;
  Variables variables_;
  std::vector<int> channels_;

  template <typename T>
  void getFirstLocationValue(const Variable & var, T & outvalue) {
    std::vector<T> array;
    getObservation(var.group(), var.variable(), array, true);
    ASSERT(array.size() > 0);
    outvalue = array[0];
  }
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMSCALEDRADIANCE_H_
