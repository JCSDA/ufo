/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDIR_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

///
/// \brief Options applying to the determination of observation error bounds as a function
/// transmittance at model top and latitude
///
class ObsErrorBoundIRParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorBoundIRParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<std::vector<float>> obserrBoundMax{"obserr_bound_max", this};

  /// Function used to set the observation bound based on Latitude (ObsErrorFactorLatRad)
  oops::RequiredParameter<Variable> obserrBoundLat{"obserr_bound_latitude", this};

  /// Function used to set the observation bound based on transmittance at model top
  /// (ObsErrorFactorTransmitTopRad)
  oops::RequiredParameter<Variable> obserrBoundTransmittop{"obserr_bound_transmittop", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};

  /// Parameter for original observation error
  oops::OptionalParameter<std::vector<float>> obserrOriginal{"error parameter vector", this};
};

///
/// \brief Determine the observation error bound (Residual Threshold) for gross check
/// as a function of transmittance at model top and latutude.
/// Errobs0   = un-inflated observation error
/// ErrobsMax = maximum observation error bound
/// Errflat   = error factor as a function of latitude
/// Errtaotop = error factor as a function of transmittance at model top
/// Residual Threshold = MIN( (3.0 * ( 1 / Errflat )^2 * (1 / Errftaotop )^2), ErrobsMax )
/// Filter out data if |obs-h(x)| > Residual Threshold
///
class ObsErrorBoundIR : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorBoundIR(const eckit::LocalConfiguration &);
  ~ObsErrorBoundIR();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorBoundIRParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDIR_H_
