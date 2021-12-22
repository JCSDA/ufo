/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SOLARZENITH_H_
#define UFO_FILTERS_OBSFUNCTIONS_SOLARZENITH_H_

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Configuration options of SolarZenith.
class SolarZenithParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SolarZenithParameters, Parameters)

 public:
  /// Set this option to `true` to skip calculations and produce missing values at locations
  /// where all simulated variables have been rejected. Default: `false`.
  oops::Parameter<bool> skipRejected{"skip rejected", false, this};
};

/// \brief Compute the solar zenith angle of observations (in degrees) as a function of their time
/// and location.
///
/// References:
/// * `Ops_Solar_Zenith` (subroutine in the Met Office OPS system): original source code
/// * Air Almanac: useful for checking GHA and DECL
/// * Norton's Star Atlas: for equation of time
/// * Robinson N., Solar Radiation, Ch. 2: for useful introduction to theory/terminology.
class SolarZenith : public ObsFunctionBase<float> {
 public:
  explicit SolarZenith(const eckit::LocalConfiguration &conf);

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override;

 private:
  SolarZenithParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SOLARZENITH_H_
