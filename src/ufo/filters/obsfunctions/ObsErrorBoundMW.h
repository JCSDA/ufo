/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDMW_H_

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
/// latitude, terrain height, and transmittance at the model top
///
class ObsErrorBoundMWParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorBoundMWParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the observation error factor applies
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// The maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<std::vector<float>> obserrBoundMax{"obserr_bound_max", this};

  /// Function to set the observation bound based on latitude
  oops::RequiredParameter<Variable> obserrBoundLat{"obserr_bound_latitude", this};

  /// Function to set the observation bound based on transmittance at model top
  oops::RequiredParameter<Variable> obserrBoundTransmittop{"obserr_bound_transmittop", this};

  /// Function to set the observation bound based on terrain height
  oops::RequiredParameter<Variable> obserrBoundTopo{"obserr_bound_topo", this};

  /// Factor applied to the derived top bound. It is 3.0 by default;
  oops::OptionalParameter<float> thresholdfactor{"threshold", this};

  /// Function to estimate observation error based on symmetric cloud amount
  // oops::RequiredParameter<Variable> obserrFunction{"obserr_function", this};
  oops::OptionalParameter<Variable> obserrFunction{"obserr_function", this};

  /// Parameter for original observation error
  oops::OptionalParameter<std::vector<float>> obserrOriginal{"error parameter vector", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
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
class ObsErrorBoundMW : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorBoundMW(const eckit::LocalConfiguration &);
  ~ObsErrorBoundMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorBoundMWParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDMW_H_
