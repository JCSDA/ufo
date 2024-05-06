/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELQUAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELQUAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsErrorModelQuad ObsFunction
class ObsErrorModelQuadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorModelQuadParameters, Parameters)

 public:
  /// x variable of the piece-wise function
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  /// channels for which to calculate the ObsError for radiances
  /// Omit channels for application to a single non-radiance filter variable.
  /// When multiple "filter variables" are provided without channels,
  /// they will have the same observation error.
  oops::OptionalParameter<std::string> chlist{"channels", this};
  /// curvature of the quadratic function
  oops::RequiredParameter<std::vector<float>> a{"a", this};
  /// x-coordinate of the quadratic function apex
  oops::RequiredParameter<std::vector<float>> b{"b", this};
  /// y-coordinate of the lower piecewise inflection point
  oops::RequiredParameter<std::vector<float>> err0{"err0", this};
  /// y-coordinate of the upper piecewise inflection point
  oops::RequiredParameter<std::vector<float>> err1{"err1", this};
  /// whether to save the xfunc values to the ObsSpace
  oops::Parameter<bool> save{"save", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Parameterize the observation error as a
///        piece-wise quadratic function of a ufo::Variable
///
/// The piece-wise function includes:
/// -# initial constant value
/// -# quadratic or inverse-quadratic growth
/// -# final constant value
///
/// and is specified by a, b, err0, and err1 as follows:
/// ~~~~
/// p0 * x^2 + p1 * x + p2 ≡ a * (x - b)^2 + c
///
/// a = p0
/// b = p1 / 2p0
/// c = p2 - p1^2 / (4 * p0)
///
/// For a < 0
/// c = err1 |-   -  _,.-----
///          |     .'   '
///          |   ,'     '
///     err0 |__/       '
///          |          '
///          '----------+-----
///                     '
///                     b
/// For a > 0
///     err1 |-   -   - ,-----
///          |         /
///          |      _.'
/// c = err0 |___.-'
///          |  .
///          '--+-------------
///             '
///             b
/// ~~~~
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - Filter: {Filter Name}
///
/// #### ABI/AHI ####
///
///       filter variables:
///       - name: brightnessTemperature
///         channels: &errassignchan 8-10
///       action:
///         name: assign error
///         error function:
///           name: ObsFunction/ObsErrorModelQuad
///           channels: *errassignchan
///           options:
///             channels: *errassignchan
///             xvar:
///               name: ObsFunction/OkamotoSCIforIR
///               channels: *errassignchan
///               options:
///                 channels: *errassignchan
///             a: [-0.069, -0.045, -0.032]
///             b: [15.0,  20.0,  25.0]
///             err0: [ 2.5,  3.2,  3.2]
///             err1: [17.0, 20.5, 21.1]
///            {save: true}
///
class ObsErrorModelQuad : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorModelQuad";}

  explicit ObsErrorModelQuad(const eckit::LocalConfiguration &);
  ~ObsErrorModelQuad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsErrorModelQuadParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELQUAD_H_
