/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_BGDDEPARTUREANOMALY_H_
#define UFO_FILTERS_OBSFUNCTIONS_BGDDEPARTUREANOMALY_H_

#include <string>
#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options for calculating the background departure anomaly between two channels.
///
class BgdDepartureAnomalyParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BgdDepartureAnomalyParameters, Parameters)

 public:
  /// \brief Lower frequency channel number
  ///
  /// Example: AMSR2 channel 11 corresponding to 36.5 GHz H-pol
  ///
  ///          channel_low_freq: 11
  oops::RequiredParameter<int> obslow{"channel_low_freq", this, {oops::minConstraint(1)}};

  /// \brief Higher frequency channel number
  ///
  /// Example: AMSR2 channel 13 corresponding to 89 GHz H-pol
  ///
  ///          channel_high_freq: 13
  oops::RequiredParameter<int> obshigh{"channel_high_freq", this, {oops::minConstraint(1)}};

  /// \brief Name of the bias correction group used to apply correction to ObsValue.
  ///
  /// Default (missing optional parameter) applies no bias
  ///
  /// Example: use ObsBias correction values
  ///
  ///          ObsBias: ObsBias
  oops::Parameter<std::string> ObsBias{"ObsBias", "", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"test_hofx", "HofX", this};
};

///
/// \brief
/// Hydrometeors scatter radiation more efficiently with increasing microwave frequency.
/// For example, at 89 GHz the scattering due to clouds is sufficient to decrease the top
///  of atmosphere brightness temperature relative to clear skies. At 36 GHz the scattering
///  is reduced, and clouds appear warm relative to surface emissions.
///
/// The anomaly (O-B minus the mean O-B) difference between these two frequencies
/// is used to diagnose cloudy scenes.
///

class BgdDepartureAnomaly : public ObsFunctionBase<float> {
 public:
  explicit BgdDepartureAnomaly(const eckit::LocalConfiguration & = eckit::LocalConfiguration());
  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  BgdDepartureAnomalyParameters options_;
  std::vector<int> channels_ = {0, 0};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_BGDDEPARTUREANOMALY_H_
