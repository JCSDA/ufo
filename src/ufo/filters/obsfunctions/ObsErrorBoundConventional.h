/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDCONVENTIONAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDCONVENTIONAL_H_

#include <memory>
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
class ObsErrorBoundConventionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorBoundConventionalParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> obsvar{"obsvar", "stationPressure", this};
  /// The maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<float> obserrBoundMax{"obserr_bound_max", this};
  /// The maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<float> obserrBoundMin{"obserr_bound_min", this};
  /// The maximum value of the observation error bound for each channel in channelList
  oops::RequiredParameter<float> obserrBoundFactor{"obserr_bound_factor", this};
  /// Function to set the observation bound based on latitude
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
class ObsErrorBoundConventional : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorBoundConventional(const eckit::LocalConfiguration &);
  ~ObsErrorBoundConventional();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
//  ObsErrorBoundConventionalParameters options_;
  std::unique_ptr<ObsErrorBoundConventionalParameters> options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORBOUNDCONVENTIONAL_H_
