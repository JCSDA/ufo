/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SATWINDINDIVERRORS_H_
#define UFO_FILTERS_OBSFUNCTIONS_SATWINDINDIVERRORS_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling SatwindIndivErrors ObsFunction
class SatwindIndivErrorsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatwindIndivErrorsParameters, Parameters)

 public:
  /// Vector error estimate addition
  oops::RequiredParameter<float> eu_add{"verror add", this};
  /// Vector error estimate multiply
  oops::RequiredParameter<float> eu_mult{"verror mult", this};
  /// Profile we are calculating error for
  oops::RequiredParameter<std::string> profile{"wind component", this};
  /// vertical coordinate to use
  oops::RequiredParameter<std::string> vcoord{"vertical coordinate", this};
  /// default pressure error (Pa)
  oops::RequiredParameter<float> default_err_p{"default pressure error", this};
  /// ignore contribution above height of minimum pressure (Pa)
  oops::OptionalParameter<float> min_press{"minimum pressure", this};
};

// -----------------------------------------------------------------------------

///
/// \brief Function calculates individual observation errors for Satwind u and v winds
///  dependent on an input height (pressure) error estimate and the wind shear.
///

class SatwindIndivErrors : public ObsFunctionBase<float> {
 public:
  explicit SatwindIndivErrors(const eckit::LocalConfiguration &);
  ~SatwindIndivErrors();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  SatwindIndivErrorsParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SATWINDINDIVERRORS_H_
