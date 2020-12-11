/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SATWINDSSPDBCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_SATWINDSSPDBCHECK_H_

#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Two options are required for this function: error_min and error_max,
///        which are bounding the ObsError variable following a GSI fix file.
///        The specification of group name for ObsError is optional but assigned
///        ObsErrorData by default.  The HofX group name is also optional to override
///        using testHofX for testing purposes (set to GsiHofX, for example).
///
class SatWindsSPDBCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatWindsSPDBCheckParameters, Parameters)

 public:
  /// the existence of min,max error values are required
  oops::RequiredParameter<float> error_min{"error_min", this};
  oops::RequiredParameter<float> error_max{"error_max", this};
  /// Name of the HofX group used to replace the default (HofX) group
  oops::Parameter<std::string> testHofX{"testHofX", "HofX", this};
  /// Name of the ObsError group used to replace the default (ObsErrorData) group
  oops::Parameter<std::string> original_obserr{"original_obserr", "ObsErrorData", this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the speed background difference (SPDB) of observation-model
///        SatWinds.  If the SPDB is less than zero, then compute a final value
///        of residual over ObsError.  That calculation is then used versus a
///        threshold value in Bounds Check filter to flag obs locations as bad.
///        In regular usage, the testHofX option would be omitted in order than HofX
///        values are used for model wind components.  Also, an optional group
///        name for the ObsError variable can be supplied if different from ObsErrorData,
///        and its value is bounded by parameters error_min and error_max.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Bounds Check
///       filter variables:
///       - name: eastward_wind
///       - name: northward_wind
///       test variables:
///       - name: SatWindsSPDBCheck@ObsFunction
///         options:
///           error_min: 1.4
///           error_max: 20.0
///       maxvalue: 1.75            # gross error * 0.7
///
class SatWindsSPDBCheck : public ObsFunctionBase {
 public:
  static const std::string classname() {return "SatWindsSPDBCheck";}

  explicit SatWindsSPDBCheck(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~SatWindsSPDBCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  SatWindsSPDBCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SATWINDSSPDBCHECK_H_
