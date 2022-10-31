/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_EMISSIVITY_DIFF_GMI_H_
#define UFO_FILTERS_OBSFUNCTIONS_EMISSIVITY_DIFF_GMI_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Compare GMI emissivity regressions over ocean.
///
class Emissivity_Diff_GMIParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(Emissivity_Diff_GMIParameters, Parameters)

 public:
  /// One of GMI channels 1-13.
  oops::OptionalParameter<int> channel{"channel", this};
  oops::OptionalParameter<float> regression_constant_1{"regression_constant_1", this};
  oops::OptionalParameter<float> regression_constant_2{"regression_constant_2", this};
  oops::OptionalParameter<std::vector<float>> regression_coeff_1{"regression_coeff_1", this};
  oops::OptionalParameter<std::vector<float>> regression_coeff_2{"regression_coeff_2", this};
};

class Emissivity_Diff_GMI : public ObsFunctionBase<float> {
 public:
  explicit Emissivity_Diff_GMI(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~Emissivity_Diff_GMI();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  Emissivity_Diff_GMIParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_EMISSIVITY_DIFF_GMI_H_
