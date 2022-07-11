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

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

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
  /// String containing the name of the wind component we are calculating the error for
  oops::RequiredParameter<std::string> profile{"wind component", this};
  /// String containing the observation vertical coordinate
  oops::RequiredParameter<std::string> obs_vcoord{"observation vertical coordinate", this};
  /// String containing the vertical coordinate to use for the model wind component
  oops::RequiredParameter<std::string> vcoord{"vertical coordinate", this};
  /// Ignore contribution above height of minimum pressure (Pa)
  oops::Parameter<float> min_press{"minimum pressure", 10000.0, this};
  /// Name of the variable containing the input height error estimates (Pa)
  oops::RequiredParameter<Variable> pressure_error{"pressure error", this};
  /// Name of the variable containing quality index values for use in the vector error calculation
  oops::RequiredParameter<Variable> quality_index{"quality index", this};
};

// -----------------------------------------------------------------------------

///
/// \brief Function calculates individual observation errors for Satwind u and v winds
///  dependent on an input pressure error estimate and the model wind shear.
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
