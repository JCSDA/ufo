/*
 * (C) Copyright 2022 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
*/

#ifndef UFO_FILTERS_OBSFUNCTIONS_SCATWINDSAMBIGUITYCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_SCATWINDSAMBIGUITYCHECK_H_

#include <string>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief An optional parameter to override the source of HofX wind components,
///        and an optional parameter for minimum wind components (default=0.5 m/s).
///
class ScatWindsAmbiguityCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ScatWindsAmbiguityCheckParameters, Parameters)

 public:
  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> test_hofx{"test_hofx", "HofX", this};
  /// minimum_uv default value set to minimum allowable value 0.0001
  oops::Parameter<float> minimum_uv{"minimum_uv", 0.0001, this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the fit of the observed wind to the model background as the length
///        of the vector difference, and compare to a hypothetical observed wind that
///        is rotated 180-degrees. Returns the difference between these values.
///        For application in ScatWinds QC, a returned positive value indicates that
///        a wind pointing in the opposite direction is a better fit to the model, which
///        implies the wrong scatterometer wind ambiguity was selected and the observation
///        should be rejected.  A threshold value of zero can be used as the max_value
///        in a Bounds Check filter to properly QC scatterometer winds.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Bounds Check
///       filter variables:
///       - name: windEastward
///       - name: windNorthward
///       test variables:
///       - name: ObsFunction/ScatWindsAmbiguityCheck
///         options:
///           test_hofx: GsiHofX
///       maxvalue: 0.
///
class ScatWindsAmbiguityCheck : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ScatWindsAmbiguityCheck";}

  explicit ScatWindsAmbiguityCheck(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~ScatWindsAmbiguityCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ScatWindsAmbiguityCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SCATWINDSAMBIGUITYCHECK_H_
