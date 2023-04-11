/*
 * (C) Copyright 2023 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
*/

#ifndef UFO_FILTERS_OBSFUNCTIONS_SATWINDSERRNORMCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_SATWINDSERRNORMCHECK_H_

#include <string>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief An optional parameter for minimum wind components
///      (default=0.0001 m/s).
class SatWindsErrnormCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatWindsErrnormCheckParameters, Parameters)

 public:
  /// minimum_uv default value set to minimum allowable value 0.0001
  oops::Parameter<float> minimum_uv{"minimum_uv", 0.0001, this};
};
// -----------------------------------------------------------------------------

/// \brief Compute and return the ratio of the expected error to the observed
///        wind speed. For application in SatWinds QC, a returned value larger
///        than a threshold value (typically 0.9) indicates that the wind
///        observation is ambiguous and should be rejected. The threshold value
///        can be used as the maxvalue in a Bounds Check filter.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Bounds Check
///       filter variables:
///       - name: windEastward
///       - name: windNorthward
///       test variables:
///       - name: ObsFunction/SatWindsErrnormCheck
///       maxvalue: 0.9
///
class SatWindsErrnormCheck : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "SatWindsErrnormCheck";}

  explicit SatWindsErrnormCheck(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~SatWindsErrnormCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  SatWindsErrnormCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SATWINDSERRNORMCHECK_H_
