/*
 * (C) Crown copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CIRCULARDIFFERENCE_H_
#define UFO_FILTERS_OBSFUNCTIONS_CIRCULARDIFFERENCE_H_

#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Parameters for the CircularDifference ObsFunction
class CircularDifferenceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CircularDifferenceParameters, Parameters)

 public:
  /// Start variable (value to subtract from in the difference calculation)
  oops::RequiredParameter<Variable> variableStart{"variable start", this};

  /// End variable (value to subtract in the difference calculation)
  oops::RequiredParameter<Variable> variableEnd{"variable end", this};

  /// Whether to return signed difference (true) or absolute difference (false)
  /// Default is true (signed difference).
  oops::Parameter<bool> computeSignedDifference{"signed", true, this};

  /// Circular period (e.g., 360.0 for degrees, 2*pi for radians, 24 for hours
  /// in a day, etc.)
  oops::RequiredParameter<float> circularPeriod{"circular period", this};
};

/// \brief Compute minimum circular difference between two variables
///
/// This ObsFunction computes the **minimum (shortest-path)** difference
/// between two values defined on a circular (cyclic) domain, accounting
/// for wrap-around at the circular period boundary.
///
/// Consider two angles: For A (variable start) and B (variable end), the
/// circular difference is the shortest angular distance from A to B.
///
/// If signed = true:
///   Returns a signed difference in the range [-period/2, period/2)
///   For angles, positive values indicate B is clockwise from A
///   For angles, negative values indicate B is counter-clockwise from A
///
/// If signed = false:
///   Returns the absolute (unsigned) difference in the range [0, period/2]
///
/// Example for wind directions (period = 360):
///   - A=10°, B=350°: signed=-20°, unsigned=20°
///   - A=350°, B=10°: signed=20°, unsigned=20°
///   - A=0°, B=180°: signed=-180°, unsigned=180°
///
/// Example for days of the week (period = 7):
///   - A=1 (Monday), B=6 (Saturday): signed=-2, unsigned=2
///   - A=3 (Wednesday), B=4 (Thursday): signed=1, unsigned=1
///
/// Usage example:
/// \code{.yaml}
///   obs function:
///     name: ObsFunction/CircularDifference
///     options:
///       variable start: ObsValue/windDirectionA
///       variable end: ObsValue/windDirectionB
///       signed: true
///       circular period: 360.0
/// \endcode
///
/// \note
///   This function always returns the minimum (shortest-path) circular
///   difference. Selecting the complementary (long-path) difference is
///   not currently supported, but may be added as a future option.
class CircularDifference : public ObsFunctionBase<float> {
 public:
  static const std::string classname() { return "CircularDifference"; }

  /// Constructor takes configuration and validates parameters
  explicit CircularDifference(const eckit::LocalConfiguration &);

  /// Destructor
  ~CircularDifference();

  /// Compute the circular difference
  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;

  /// Return the required variables for this ObsFunction
  const ufo::Variables &requiredVariables() const;

 private:
  /// Parameters for this ObsFunction
  CircularDifferenceParameters options_;

  /// Required variables
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CIRCULARDIFFERENCE_H_
