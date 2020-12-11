/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELSTEPWISELINEAR_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELSTEPWISELINEAR_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsErrorModelStepwiseLinear ObsFunction
class ObsErrorModelStepwiseLinearParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorModelStepwiseLinearParameters, Parameters)

 public:
  /// the name of the xvar
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  /// vector of X-steps
  oops::RequiredParameter<std::vector<float>> xvals{"xvals", this};
  /// vector of error values corresponding to vector of X-steps
  oops::RequiredParameter<std::vector<float>> errors{"errors", this};
};

// -----------------------------------------------------------------------------

/// \brief Parameterize the observation error as a series of steps with
///        linear interpolation between each step.
///
/// This routine was designed to mimic the GSI fix-file of prepobs_errtable.txt
/// Input is a vector of x-values (e.g. pressures) and corresponding vector of obserrors.
/// Interpolation in X-coordinate requires the value of X for which the output, Y, is
/// calculated using linear interp of obserrors between the steps.
/// ~~~~
///
///          + err_n                             o
///          |                                /
///          |                             /
///          + err_n-1                  o
///          |                       /
///  obserr  +                     *
///          |                  /
///          |               /
///     err2 +            o
///          |          /
///          |       /
///     err1 +     o
///          |  /
///     err0 +-----+------+---------*---+--------+
///          p_0, p_1,   p_2,...   p*  p_n-1,   p_n
///
/// ~~~~
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - Filter: {Filter Name}
///
/// #### Example for air temperature assigned obserror by pressure (e.g., sonde data) ####
///
///       filter variables:
///       - name: air_temperature
///       action:
///         name: assign error
///         error function:
///           name: ObsErrorModelStepwiseLinear@ObsFunction
///           options:
///             xvar:
///               name: air_pressure@ObsValue
///             xvals: [110000, 85000, 50000, 25000, 10000, 1]   #Pressure (Pa)
///             errors: [1.1, 1.3, 1.8, 2.4, 4.0, 4.5]
///
class ObsErrorModelStepwiseLinear : public ObsFunctionBase {
 public:
  static const std::string classname() {return "ObsErrorModelStepwiseLinear";}

  explicit ObsErrorModelStepwiseLinear(const eckit::LocalConfiguration);
  ~ObsErrorModelStepwiseLinear();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsErrorModelStepwiseLinearParameters options_;
  bool isAscending_ = true;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELSTEPWISELINEAR_H_
