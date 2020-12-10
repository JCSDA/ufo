/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SATWINDSLNVDCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_SATWINDSLNVDCHECK_H_

#include <string>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Only one option to override the source of HofX wind components.
///
class SatWindsLNVDCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatWindsLNVDCheckParameters, Parameters)

 public:
  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"testHofX", "HofX", this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the log-normal vector difference (LNVD) of observation-model
///        SatWinds.  If the LNVD is greater than a threshold (usually 3.0), then
///        flag the location as bad.  In regular usage, the testHofX option would
///        be omitted in order that HofX is used for the model wind components.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Bounds Check
///       filter variables:
///       - name: eastward_wind
///       - name: northward_wind
///       test variables:
///       - name: SatWindsLNVDCheck@ObsFunction
///         options:
///           testHofX: GsiHofX
///       maxvalue: 3
///
class SatWindsLNVDCheck : public ObsFunctionBase {
 public:
  static const std::string classname() {return "SatWindsLNVDCheck";}

  explicit SatWindsLNVDCheck(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~SatWindsLNVDCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  SatWindsLNVDCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SATWINDSLNVDCHECK_H_
