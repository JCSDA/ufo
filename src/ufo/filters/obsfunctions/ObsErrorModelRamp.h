/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELRAMP_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELRAMP_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Options controlling the ObsErrorModelRamp ObsFunction
class ObsErrorModelRampParameters : public oops::Parameters {
 public:
  /// x variable of the piece-wise function
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  /// channels; channels for which to calculate the ObsError for radiances
  /// Omit channels for application to a single non-radiance filter variable.
  /// When multiple "filter variables" are provided without channels,
  /// they will have the same observation error.
  oops::OptionalParameter<std::string> chlist{"channels", this};
  /// x-coordinate of the lower ramp inflection point
  oops::RequiredParameter<std::vector<float>> x0{"x0", this};
  /// x-coordinate of the upper ramp inflection point
  oops::RequiredParameter<std::vector<float>> x1{"x1", this};
  /// y-coordinate of the lower ramp inflection point
  oops::RequiredParameter<std::vector<float>> err0{"err0", this};
  /// y-coordinate of the upper ramp inflection point
  oops::RequiredParameter<std::vector<float>> err1{"err1", this};
  /// whether to save xvar values to the ObsSpace
  oops::Parameter<bool> save{"save", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Parameterize the observation error as a
///        piece-wise linear function of a ufo::Variable
///
/// The output function is specified by the coordinates of the
/// inflection points and includes:
/// -# initial constant value
/// -# linear growth or decay (ramp)
/// -# final constant value
///
/// ~~~~
/// err1 |-   -   - *---
///      |        ,'.
///      |      ,'  .
///      |    ,'    .
/// err0 |--*'      .
///      |  .       .
///      '--+-------+---
///         '       '
///        x0      x1
/// ~~~~
///
/// Notes:
/// - for a decaying ramp, set err1 < err0
/// - for a step function, set x0 == x1
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - Filter: {Filter Name}
///
/// #### AMSUA ####
///
///       filter variables:
///       - name: brightness_temperature
///         channels: &errassignchan 1-15
///       action:
///         name: assign error
///         error function:
///           name: ObsErrorModelRamp@ObsFunction
///           channels: *errassignchan
///           options:
///             channels: *errassignchan
///             xvar:
///               name: CLWRetMean@ObsFunction
///               options:
///                 clwret_ch238: 1
///                 clwret_ch314: 2
///                 clwret_types: [ObsValue, HofX]
///                 bias_application: HofX
///             x0:    [ 0.050,  0.030,  0.030,  0.020,  0.000,
///                      0.100,  0.000,  0.000,  0.000,  0.000,
///                      0.000,  0.000,  0.000,  0.000,  0.030]
///             x1:    [ 0.600,  0.450,  0.400,  0.450,  1.000,
///                      1.500,  0.000,  0.000,  0.000,  0.000,
///                      0.000,  0.000,  0.000,  0.000,  0.200]
///             err0: [ 2.500,  2.200,  2.000,  0.550,  0.300,
///                     0.230,  0.230,  0.250,  0.250,  0.350,
///                     0.400,  0.550,  0.800,  3.000,  3.500]
///             err1: [20.000, 18.000, 12.000,  3.000,  0.500,
///                     0.300,  0.230,  0.250,  0.250,  0.350,
///                     0.400,  0.550,  0.800,  3.000, 18.000]
///            {save: true}
///
/// #### ABI/AHI ####
///
///       filter variables:
///       - name: brightness_temperature
///         channels: &errassignchan 8-10
///       action:
///         name: assign error
///         error function:
///           name: ObsErrorModelRamp@ObsFunction
///           channels: *errassignchan
///           options:
///             channels: *errassignchan
///             xvar:
///               name: SymmCldImpactIR@ObsFunction
///               channels: *errassignchan
///               options:
///                 channels: *errassignchan
///             x0: [ 0.0,  0.0,  1.0]
///             x1: [15.0, 20.0, 25.0]
///             err0: [ 2.5,  3.2,  3.2]
///             err1: [17.0, 20.5, 21.1]
///
/// #### Non-radiance ObsTypes ####
///
///       filter variables:
///       - name: {filter variable name}
///       action:
///         name: assign error
///         error function:
///           name: ObsErrorModelRamp@ObsFunction
///           options:
///             xvar:
///               name: {xvar}@[ObsFunction, GeoVaLs, ObsDiag, ObsValue, etc...]
///               options: {xvar options}
///             x0: [{X0}]
///             x1: [{X1}]
///             err0: [{ERR0}]
///             err1: [{ERR1}]
///
class ObsErrorModelRamp : public ObsFunctionBase {
 public:
  static const std::string classname() {return "ObsErrorModelRamp";}

  explicit ObsErrorModelRamp(const eckit::LocalConfiguration);
  ~ObsErrorModelRamp();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsErrorModelRampParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELRAMP_H_
