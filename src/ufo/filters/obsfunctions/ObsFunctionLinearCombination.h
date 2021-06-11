/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsFunctionLinearCombination ObsFunction
class LinearCombinationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearCombinationParameters, Parameters)

 public:
  /// Input variables of the linear combination
  oops::RequiredParameter<std::vector<ufo::Variable>> variables{"variables", this};
  /// coefficient associated with the above variables
  oops::RequiredParameter<std::vector<float>> coefs{"coefs", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs a linear combination of variables
///
/// For example, the following
///
///  obs function:
///    name: LinearCombination@ObsFunction
///    options:
///      variables: [representation_error@GeoVaLs,
///                  sea_water_temperature@ObsError]
///      coefs: [0.1,
///              1.0]
///
/// will return 0.1 *representation_error@GeoVaLs +
///             1.0 * sea_water_temperature@ObsError
///


class LinearCombination : public ObsFunctionBase {
 public:
  explicit LinearCombination(const eckit::LocalConfiguration &);
  ~LinearCombination();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  LinearCombinationParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_
