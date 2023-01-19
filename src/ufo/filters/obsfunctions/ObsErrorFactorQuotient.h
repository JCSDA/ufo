/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORQUOTIENT_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORQUOTIENT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorQuotient ObsFunction
class ObsErrorFactorQuotientParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorQuotientParameters, Parameters)

 public:
  /// the name of the numerator and denominator variables (with group names).
  oops::RequiredParameter<Variable> numerator{"numerator", this};
  oops::RequiredParameter<Variable> denominator{"denominator", this};
  oops::Parameter<bool> save{"save", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Calculate the ratio of two variables, typically related to ObsError
///
/// This routine was designed to mimic the GSI-Observer method of rejecting obs data
/// when the ratio of the final ObsError, after inflation, is more than a threshold
/// amount greater than the starting ObsError.  A boolean optional save variable
/// is available, false by default.  The ObsFunction is simple division expected to
/// be used with Bounds Check above a maximum threshold.
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - Filter: Bounds Check
///       filter variables:
///       - name: airTemperature
///       action:
///         name: reject
///       maxvalue: 3.6
///       test variables:
///       - name: ObsFunction/ObsErrorFactorQuotient
///           options:
///             numerator:
///               name: ObsErrorData/airTemperature   # After inflation step
///             denominator:
///               name: ObsError/airTemperature
///       defer to post: true                         # Likely necessary for order of filters
///
class ObsErrorFactorQuotient : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorFactorQuotient";}

  explicit ObsErrorFactorQuotient(const eckit::LocalConfiguration);
  ~ObsErrorFactorQuotient();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsErrorFactorQuotientParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORQUOTIENT_H_
