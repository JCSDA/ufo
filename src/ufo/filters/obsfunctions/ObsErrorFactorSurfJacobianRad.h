/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSURFJACOBIANRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSURFJACOBIANRAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to observation error inflation as a function
//  of surface temperature jacobian and surface emissivity jacobian
///
class ObsErrorFactorSurfJacobianRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorSurfJacobianRadParameters, Parameters)

 public:
  /// List of channels to which the observation error factor applies
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the observation error factor applies
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// Observation error scale factors applied to surface temperature jacobians
  /// over five surface types: [sea, land, ice, snow and mixed]
  oops::RequiredParameter<std::vector<float>> obserrScaleFactorTsfc{"obserr_dtempf", this};

  /// Observation error scale factors applied to surface emissivity jacobians
  /// over five surface types: [sea, land, ice, snow and mixed]
  oops::RequiredParameter<std::vector<float>> obserrScaleFactorEsfc{"obserr_demisf", this};

  /// Bias term usage indicator
  oops::OptionalParameter<bool> useBiasTerm{"use_biasterm", this};

  /// Name of the group for bias correction terms used to replace the default group
  /// (default is ObsBiasTerm)
  oops::Parameter<std::string> testBiasTerm{"test_biasterm", "ObsBiasTerm", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

///
/// \brief Error Inflation Factor (EIF) as a function of weighted surface temperature
/// jacobian and surface emissivity jacobian
/// Jtemp  = surface temperature jacobian
/// Jemis  = surface emissivity jacobian
/// Wtemp  = empirical constant as a function of surface type applied to Jtemp
/// Wemis  = empirical constant as a function of surface type applied to Jemis
/// Beta   = ( Wtemp * Jtemp + Wemis * Jemis )^2
/// Errinv = inverse of effective observation error variance
/// EIF    = SQRT [ 1 / ( 1 / (1 + Errinv * Beta) ]
///
class ObsErrorFactorSurfJacobianRad : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorSurfJacobianRad(const eckit::LocalConfiguration &);
  ~ObsErrorFactorSurfJacobianRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorFactorSurfJacobianRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSURFJACOBIANRAD_H_
