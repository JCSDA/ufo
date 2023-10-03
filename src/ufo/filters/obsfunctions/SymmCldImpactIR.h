/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACTIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACTIR_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Options controlling Symmetric Cloud Impact for IR instruments
class SymmCldImpactIRParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SymmCldImpactIRParameters, Parameters)

 public:
  /// channels for which SCI will be calculated
  oops::RequiredParameter<std::string> chlist{"channels", this};
  oops::Parameter<bool> scale_by_omb{"scale by omb", false, this};
  oops::Parameter<float> sigmoid_c1{"sigmoid constant 1", 10.0f, this};
  oops::Parameter<float> sigmoid_c2{"sigmoid constant 2", 10.0f, this};
};

// -----------------------------------------------------------------------------

/// \brief Okamoto et al. symmetric cloud impact (SCI) function
///
/// Okamoto, K., McNally, A.P. and Bell, W. (2014), Progress towards the
///   assimilation of all‐sky infrared radiances: an evaluation of cloud
///   effects. Q.J.R. Meteorol. Soc., 140: 1603-1614. doi:10.1002/qj.2242
class SymmCldImpactIR : public ObsFunctionBase<float> {
 public:
  explicit SymmCldImpactIR(const eckit::LocalConfiguration);
  ~SymmCldImpactIR();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  SymmCldImpactIRParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SYMMCLDIMPACTIR_H_
