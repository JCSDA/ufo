/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDRAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

///
/// \brief Options applying to the determination of observation error bounds as a function
/// transmittance at model top and latitude
///
class ObsErrorBoundRadParameters : public oops::Parameters {
 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<std::vector<float>> obserrBoundMax{"obserr_bound_max", this};

  /// Name of the variable used to retrieve the cloud liquid water from observation
  oops::RequiredParameter<Variable> obserrBoundLat{"obserr_bound_latitude", this};

  /// Name of the variable used to retrieve the cloud liquid water from background
  oops::RequiredParameter<Variable> obserrBoundTransmittop{"obserr_bound_transmittop", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"obserr_test", "ObsErrorData", this};
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
class ObsErrorBoundRad : public ObsFunctionBase {
 public:
  explicit ObsErrorBoundRad(const eckit::LocalConfiguration &);
  ~ObsErrorBoundRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorBoundRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDRAD_H_
