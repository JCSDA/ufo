/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSITUDEPENDMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSITUDEPENDMW_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

///
/// \brief Options applying to the situation-dependent error inflation factor
///
class ObsErrorFactorSituDependMWParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorSituDependMWParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the observation error factor applies
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// Observation error for each channel under the clear-sky condition
  oops::RequiredParameter<std::vector<float>> obserrClearSky{"obserr_clearsky", this};

  /// Function to retrieve the cloud liquid water from the observation
  oops::RequiredParameter<Variable> clwobsFunction{"clwobs_function", this};

  /// Function to retrieve the cloud liquid water from the simulated observation
  oops::RequiredParameter<Variable> clwbkgFunction{"clwbkg_function", this};

  /// Function to retrieve the scattering index from the observation
  oops::RequiredParameter<Variable> scatobsFunction{"scatobs_function", this};

  /// Function to get the cloud match index based on cloud amount retrieved from
  /// background and observation
  oops::RequiredParameter<Variable> clwmatchidxFunction{"clwmatchidx_function", this};

  /// Function used to estimate observation error based on symmetric cloud amount
  oops::OptionalParameter<Variable> obserrFunction{"obserr_function", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"test_hofx", "HofX", this};
};

///
/// \brief Situation-dependent error inflation factor based on
/// retrieved cloud liquid water from background and observation, scattering index,
/// surface wind speed, and cloud information match index over the ocean
///
class ObsErrorFactorSituDependMW : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorSituDependMW(const eckit::LocalConfiguration &);
  ~ObsErrorFactorSituDependMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  ObsErrorFactorSituDependMWParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSITUDEPENDMW_H_
