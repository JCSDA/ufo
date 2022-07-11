/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_SSMIS_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_SSMIS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Option to override varGroup default of ObsValue for all channels with ObsBias
///
class CLWRetMW_SSMISParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CLWRetMW_SSMISParameters, Parameters)

 public:
  /// Names of the data group used to retrieve the cloud liquid water
  /// Example: get retrieved channels from observation, ObsValue (default) or bias-corrected obs.
  oops::Parameter<std::string> varGroup{"varGroup", "ObsValue", this};
  oops::RequiredParameter<int> ch19h{"ch19h", this};
  oops::RequiredParameter<int> ch19v{"ch19v", this};
  oops::RequiredParameter<int> ch22v{"ch22v", this};
  oops::RequiredParameter<int> ch37h{"ch37h", this};
  oops::RequiredParameter<int> ch37v{"ch37v", this};
  oops::RequiredParameter<int> ch91v{"ch91v", this};
  oops::RequiredParameter<int> ch91h{"ch91h", this};
};

///
/// \brief Retrieve cloud liquid water from channels 12-18 of SSMIS data.
///
/// Reference:
/// Weng, F., R. R. Ferraro, and N. C. Grody,2000: "Effects of AMSU cross-scan Symmetry of
///          brightness temperatures on  retrieval of atmospheric and surface parameters",
///          Ed. P. Pampaloni and S. Paloscia, VSP, Netherlands, 255-262, 2000.
/// Yan B. and F. Weng, 'Intercalibration between Special Sensor Microwave Imager and
///          Sounder (SSMIS) and Special Sensor Microwave Imager (SSM/I)', TGARS Special
///          Issue on the DMSP SSMIS, 46, 984-995.
///
class CLWRetMW_SSMIS : public ObsFunctionBase<float> {
 public:
  explicit CLWRetMW_SSMIS(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~CLWRetMW_SSMIS();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  static void cloudLiquidWater(const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               const std::vector<float> &,
                               std::vector<float> &,
                               std::vector<float> &);

 private:
  ufo::Variables invars_;
  CLWRetMW_SSMISParameters options_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLWRETMW_SSMIS_H_
