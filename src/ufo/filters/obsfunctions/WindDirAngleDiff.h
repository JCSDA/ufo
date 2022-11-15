/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_WINDDIRANGLEDIFF_H_
#define UFO_FILTERS_OBSFUNCTIONS_WINDDIRANGLEDIFF_H_

#include <string>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief An optional parameter to override the source of HofX wind components,
///        and an optional parameter for minimum wind components (default=0.5 m/s).
///
class WindDirAngleDiffParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WindDirAngleDiffParameters, Parameters)

 public:
  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> test_hofx{"test_hofx", "HofX", this};
  oops::Parameter<float> minimum_uv{"minimum_uv", 0.5, this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the wind direction angle difference between observation and model.
///        For application in SatWinds QC, a difference larger than 50 degrees is
///        typically used to reject data.  This threshold can be used as the max_value
///        in a Bounds Check filter.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Bounds Check
///       filter variables:
///       - name: windEastward
///       - name: windNorthward
///       test variables:
///       - name: ObsFunction/WindDirAngleDiff
///         options:
///           test_hofx: GsiHofX
///       maxvalue: 50
///
class WindDirAngleDiff : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "WindDirAngleDiff";}

  explicit WindDirAngleDiff(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~WindDirAngleDiff();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  WindDirAngleDiffParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_WINDDIRANGLEDIFF_H_
