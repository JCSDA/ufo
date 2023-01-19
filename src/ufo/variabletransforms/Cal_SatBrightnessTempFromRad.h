/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_SATBRIGHTNESSTEMPFROMRAD_H_
#define UFO_VARIABLETRANSFORMS_CAL_SATBRIGHTNESSTEMPFROMRAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

enum class RadianceUnits {
  WAVENUMBER, WAVELENGTH, FREQUENCY
};

struct RadianceUnitsParameterTraitsHelper {
  typedef RadianceUnits EnumType;
  static constexpr char enumTypeName[] = "RadianceUnits";
  static constexpr util::NamedEnumerator<RadianceUnits> namedValues[] = {
    { RadianceUnits::WAVENUMBER, "wavenumber" },
    { RadianceUnits::WAVELENGTH, "wavelength" },
    { RadianceUnits::FREQUENCY, "frequency"}
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::RadianceUnits> :
    public EnumParameterTraits<ufo::RadianceUnitsParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration parameters for the radiance to brightness temperature variable transformation
class Cal_SatBrightnessTempFromRadParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_SatBrightnessTempFromRadParameters, VariableTransformParametersBase);

 public:
  /// The variable that will be transformed.  This will take the same format as other variables e.g
  /// transform variable:
  ///   name: radiance@ObsValue
  ///   channels: 1-20
  oops::RequiredParameter<Variable> transformVariable{"transform from", this};

  /// This variable contains the spectral information for each channel.  The term spectral refers to
  /// either frequency, wavelength or wavenumber for a particular channel.  This will take the same
  /// format as other variables e.g
  /// spectral variable:
  ///   name: Metadata/sensorCentralWavenumber
  ///   channels: 1-20
  oops::RequiredParameter<Variable> spectralVariable{"spectral variable", this};

  /// This variable defines the units that are used for the radiance which is technically
  /// a spectral radiance as it is a function of either frequency, wavelength or wavenumber.
  /// This has the following available options: "frequency", "wavenumber" and "wavelength"
  /// with units of Hz, m^-1 and microns, respectively.  The correct units much match
  /// with the spectral variable to have a meaningful conversion to brightness temperature.
  oops::RequiredParameter<RadianceUnits> radianceUnits{"radiance units", this};

  /// The minimum allowed value for the output brightness temperature
  oops::OptionalParameter<float> minvalue{"minimum value", this};

  /// The maximum allowed value for the output brightness temperature
  oops::OptionalParameter<float> maxvalue{"maximum value", this};

  /// The planck1 constant which equals (2*h*c*c) - (1.191042972e-16 W / (m^2.sr.m-4)).
  /// This argument is to allow for rounding differences in porting.
  oops::Parameter<double> planck1{"planck1", 1.191042972e-16, this};

  /// The planck2 constant which equals (h*c / T_b)  - (1.4387769e-2 m.K).
  /// This argument is to allow for rounding differences in porting.
  oops::Parameter<double> planck2{"planck2", 1.4387769e-2, this};
};

/*!
* \brief Convert radiance to brightness temperature.
*
* Performs a conversion of radiance to brightness temperature (K). The
* radiance can be in wavenumber (m-1), wavelength (microns) or frequency (Hz) format.
* The radiance units should match the "spectral variable" units or you
* will get incorrect results.
* The output is written to "DerivedObsValue/brightness_temperature"
*
* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variable Transforms
*   Transform: SatBrightnessTempFromRad
*   transform from:
*     name: ObsValue/radiance
*     channels: *all_channels
*   spectral variable:
*     name: MetaData/sensorCentralWavenumber
*     channels: *all_channels
*   radiance units: wavenumber
*   minimum value: 150
*   maximum value: 350
* \endcode
*
* See Cal_SatBrightnessTempFromRadParameters for filter setup.
*/

class Cal_SatBrightnessTempFromRad : public TransformBase {
 public:
  typedef Cal_SatBrightnessTempFromRadParameters Parameters_;

  Cal_SatBrightnessTempFromRad(const Parameters_ &options,
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
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_SATBRIGHTNESSTEMPFROMRAD_H_
