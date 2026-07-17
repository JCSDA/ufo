/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PERCENTILEPARAMETERS_H_
#define UFO_FILTERS_PERCENTILEPARAMETERS_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"

namespace ufo {

/// \brief Options controlling the operation of the percentile filter.
///
/// Inherits from TrackCheckUtilsParameters to gain access to
/// station_id_variable parameter which can control the ObsAccessor grouping
/// behaviour if obs grouping has not been specified in the ObsSpace.
class PercentileParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(PercentileParameters, TrackCheckUtilsParameters)

 public:
  /// Lower percentile threshold(s) for each filter variable.
  /// Data below this percentile will be rejected.
  /// Must be in range [0, 100].
  oops::OptionalParameter<std::vector<float>> lowerPercentiles{
      "lower percentiles", this};

  /// Upper percentile threshold(s) for each filter variable.
  /// Data above this percentile will be rejected.
  /// Must be in range [0, 100].
  oops::OptionalParameter<std::vector<float>> upperPercentiles{
      "upper percentiles", this};

  /// Whether the central range comparison is inclusive (<=, >=) or exclusive
  /// (<, >). For each filter variable, when true:
  ///   lower_percentile <= passing_data <= upper_percentile
  /// When false:
  ///   lower_percentile < passing_data < upper_percentile
  /// Default is true (inclusive).
  /// This can be a single value applied to all variables `[true]` or a list of
  /// values with one value per variable `[true, false, ...]`.
  oops::Parameter<std::vector<bool>> inclusiveCentralRange{
      "inclusive central range", {true}, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_PERCENTILEPARAMETERS_H_
