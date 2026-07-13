/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACT_H_
#define UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Options controlling Symmetric Cloud Impact computation
class SymmCldImpactParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SymmCldImpactParameters, Parameters)

 public:
  /// channels for which SCI will be calculated
  oops::RequiredParameter<std::string> chlist{"channels", this};

  /// Per-channel brightness temperature limit thresholds (Harnish et al. 2016).
  /// If provided, Harnish2016 formulation is used; otherwise Okamoto2014.
  oops::OptionalParameter<std::vector<float>> btlim{"btlim", this};

  /// Scale by obs-minus-background sigmoid function (Okamoto2014 only)
  oops::Parameter<bool> scale_by_omb{"scale by omb", false, this};
  /// Slope of sigmoid function (Okamoto2014 only)
  oops::Parameter<float> sigmoid_c1{"sigmoid constant 1", 10.0f, this};
  /// 50-percent value of sigmoid function (Okamoto2014 only)
  oops::Parameter<float> sigmoid_c2{"sigmoid constant 2", 10.0f, this};
};

// -----------------------------------------------------------------------------

/// \brief Symmetric cloud impact (SCI) obsfunction
///
/// Supports two formulations selected automatically based on options:
///
/// If \p btlim is provided, uses the Harnish et al. (2016) formulation:
///   Harnisch, F., M. Weissmann, and A. Perianez, 2016: Error model for the
///   assimilation of cloud-affected infrared satellite observations in an
///   ensemble data assimilation system. Q.J.R. Meteorol. Soc., 142, 1797-1808.
///   doi:10.1002/qj.2776
///
/// Otherwise, uses the Okamoto et al. (2014) formulation:
///   Okamoto, K., McNally, A.P. and Bell, W. (2014), Progress towards the
///   assimilation of all-sky infrared radiances: an evaluation of cloud
///   effects. Q.J.R. Meteorol. Soc., 140: 1603-1614. doi:10.1002/qj.2242
///
class SymmCldImpact : public ObsFunctionBase<float> {
 public:
  explicit SymmCldImpact(const eckit::LocalConfiguration);
  ~SymmCldImpact();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  SymmCldImpactParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACT_H_
