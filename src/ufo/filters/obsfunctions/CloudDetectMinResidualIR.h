/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLOUDDETECTMINRESIDUALIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLOUDDETECTMINRESIDUALIR_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to Cloud Detection Algorithm (Minimum Residual Method)
/// for Infrared sensors
///
class CloudDetectMinResidualIRParameters : public oops::Parameters {
 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Useflag (-1: not used; 0: monitoring; 1: used) for each channel in channelList
  oops::RequiredParameter<std::vector<int>> useflagChannel{"use_flag", this};

  /// Useflag (-1: not used; 1: used) indicating if the channel is used for cloud detection
  oops::RequiredParameter<std::vector<int>> useflagCloudDetect{"use_flag_clddet", this};

  /// Observation error scale factors applied to surface temperature jacobians
  /// over 5 surface types: [sea, land, ice, snow and mixed]
  oops::RequiredParameter<std::vector<float>> obserrScaleFactorTsfc{"obserr_dtempf", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"obserr_test", "ObsErrorData", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"hofx_test", "HofX", this};

  /// Name of the bias correction group used to replace the default group (default is ObsBias)
  oops::Parameter<std::string> testBias{"bias_test", "ObsBias", this};
};

///
/// \brief Cloud Detection Algorithm (Minimum Residual Method) for Infrared sensors
/// using selected channels from 15 microns CO2 absorption band
/// Output of this function:
/// 0 = channel is not affected by clouds (clear channel)
/// 1 = channel is affected by clouds (cloudy channel)
/// 2 = channel is not affected by clouds but too sensitive to surface condition
///
class CloudDetectMinResidualIR : public ObsFunctionBase {
 public:
  explicit CloudDetectMinResidualIR(const eckit::LocalConfiguration &);
  ~CloudDetectMinResidualIR();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  CloudDetectMinResidualIRParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLOUDDETECTMINRESIDUALIR_H_
