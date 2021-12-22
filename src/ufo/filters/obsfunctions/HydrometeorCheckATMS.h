/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKATMS_H_
#define UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKATMS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options controling the cloud and precipitation checks for WM sensors
///
class HydrometeorCheckATMSParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HydrometeorCheckATMSParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Observation error for each channel under clear-sky condition
  oops::RequiredParameter<std::vector<float>> obserrClearSky{"obserr_clearsky", this};

  /// Function used to estimate observation error based on symmetric cloud amount
  /// (ObsErrorModelRamp)
  oops::RequiredParameter<Variable> obserrFunction{"obserr_function", this};

  /// Function used to retrieve the cloud liquid water from observation (CLWRetMW)
  oops::RequiredParameter<Variable> clwretFunction{"clwret_function", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"test_hofx", "HofX", this};

  /// Name of the group for bias correction used to replace the default group (default is
  /// ObsBiasData)
  /// Example: use observation bias correction values from GSI
  ///          test_bias: GsiObsBias
  oops::Parameter<std::string> testBias{"test_bias", "ObsBiasData", this};

  /// Name of the group for bias correction terms used to replace the default group
  /// (default is ObsBiasTerm)
  /// Example: use observation bias correction terms from GSI
  ///          test_biasterm: GsiObsBiasTerm
  oops::Parameter<std::string> testBiasTerm{"test_biasterm", "ObsBiasTerm", this};
};

///
/// \brief Cloud and precipitation checks for ATMS
/// Checks for all observations:
///   (1) Sanity check on observaton values
///   (2) Sanity check on retrieved CLW values
/// Checks for observation over ocean include:
///   (1) Scattering check based on 54.4GHz channel
///   (2) Scattering check based on 53.6GHz channel
///   (3) Sensitivity to surface emissivity
/// Checks observation over non-ocean surface include:
///   (1) Scattering check based on 54.4GHz channel
///   (2) Thick cloud check based on 52.8GHz channel
///   (3) Sensitivity to surface emissivity
/// Output of this function:
///    0 = channel is not affected by thick clouds and precipitation
///    1 = channel is affected by thick clouds and precipitataion
///
class HydrometeorCheckATMS : public ObsFunctionBase<float> {
 public:
  explicit HydrometeorCheckATMS(const eckit::LocalConfiguration &);
  ~HydrometeorCheckATMS();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  HydrometeorCheckATMSParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_HYDROMETEORCHECKATMS_H_
