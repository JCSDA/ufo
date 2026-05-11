/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_PROFILEVERTICALSMOOTHING_H_
#define UFO_FILTERS_OBSFUNCTIONS_PROFILEVERTICALSMOOTHING_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
///
/// \brief Options: Specify the name of the variables to apply the function to,
/// the filter widths and heights to use for the smoothing, and the order of
/// the polynomial to use for the local regression.
///
class ProfileVerticalSmoothingParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProfileVerticalSmoothingParameters, Parameters)

 public:
  /// Name of variable to apply the function to
  oops::RequiredParameter<Variable> varApply{
    "variable",
    "The variable to apply the smoothing to",
    this};

  oops::RequiredParameter<Variable> heightVariable{
    "heightVariable",
    "The variable representing the height coordinate",
    this};

  oops::OptionalParameter<std::string> filterWidth{
    "filterWidth",
    "List of filter widths to apply - must be in ascending order of height."
    " Default values are 100,100.",
    this};

  oops::OptionalParameter<std::string> filterHeight{
    "filterHeight",
    "List of filter heights to apply - must be in ascending order of height."
    " Default values are 0,60000.",
    this};

  oops::Parameter<int> polynomialOrder{
    "polynomialOrder",
    "Order of polynomial to use for local regression (default is 3, i.e. cubic) "
    "Must be positive",
    3,
    this};
};

///
/// \brief Function calculates a smoothed version of the input variable by
/// applying a local polynomial regression
///
class ProfileVerticalSmoothing : public ObsFunctionBase<float> {
 public:
  explicit ProfileVerticalSmoothing(const eckit::LocalConfiguration &);
  ~ProfileVerticalSmoothing();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ProfileVerticalSmoothingParameters parameters_;
  std::vector<float> filterWidths_;
  std::vector<float> filterHeights_;
  int polyOrder_;
  // Functions used in this class
  float getFilterWidth(const std::vector<float>&, const std::vector<float>&,
                       const float) const;
  std::vector<float> parseCommaSeparatedFloats(const std::string &) const;
};
// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_PROFILEVERTICALSMOOTHING_H_
