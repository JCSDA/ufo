/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_ASSIGNVALUEEQUALCHANNELS_H_
#define UFO_FILTERS_OBSFUNCTIONS_ASSIGNVALUEEQUALCHANNELS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
///
/// \brief Options: Input variable, test value, and values to assign.
///
class AssignValueEqualChannelsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AssignValueEqualChannelsParameters, Parameters)

 public:
  /// Variable name to grab the input values from
  oops::RequiredParameter<Variable> variable{"variable", this};

  /// Value to look for in input, which will yield "good" output
  oops::Parameter<int> testValue{"testValue", 0, this};

  /// Value to give if input is "good"
  oops::Parameter<float> assignEqual{"assignEqual", 0.1, this};

  /// Value to give if input is "bad"
  oops::Parameter<float> assignNotEqual{"assignNotEqual", 1.0, this};
};

///
/// \brief Function to assign a value if the input matches a test value. If not,
///        then assign an alternative value.
///
class AssignValueEqualChannels : public ObsFunctionBase<float> {
 public:
  explicit AssignValueEqualChannels(const eckit::LocalConfiguration &);
  ~AssignValueEqualChannels();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  AssignValueEqualChannelsParameters options_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_ASSIGNVALUEEQUALCHANNELS_H_
