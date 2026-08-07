/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDSPECIFICHUMIDITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDSPECIFICHUMIDITY_H_

#include <string>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ModelHeightAdjustedSpecificHumidity ObsFunction
class ModelHeightAdjustedSpecificHumidityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightAdjustedSpecificHumidityParameters, Parameters)

 public:
  /// Observed specific humidity variable (default: ObsValue/specificHumidityAt2M)
  oops::Parameter<Variable> observedSpecificHumidity{
      "observed specific humidity variable",
      Variable("ObsValue/specificHumidityAt2M"), this};

  /// Observed temperature variable at station height (for computing RH from q_obs)
  oops::Parameter<Variable> observedTemperature{
      "observed temperature variable",
      Variable("ObsValue/airTemperatureAt2M"), this};

  /// Adjusted temperature variable (from prior Variable Assignment)
  oops::Parameter<Variable> adjustedTemperature{
      "adjusted temperature variable",
      Variable("DerivedObsValue/airTemperatureAt2M"), this};

  /// Observed station pressure variable (for computing q_sat at station height)
  oops::Parameter<Variable> stationPressure{
      "station pressure variable",
      Variable("ObsValue/stationPressure"), this};

  /// Model surface pressure variable (for computing q_sat at model surface)
  oops::Parameter<Variable> modelSurfacePressure{
      "model surface pressure variable",
      Variable("GeoVaLs/air_pressure_at_surface"), this};
};

/// \brief Adjusts specific humidity to preserve relative humidity when temperature is
/// height-adjusted from station elevation to model height.
///
/// The algorithm:
///   1. Compute RH at station: RH = q_obs / q_sat(T_obs, P_station)
///   2. Compute adjusted q at model surface: q_adj = RH * q_sat(T_adj, P_model_surface)
///
/// Where T_adj is the height-adjusted temperature (from a prior Variable Assignment using
/// ModelHeightAdjustedAirTemperature). If station pressure is missing for a location,
/// model surface pressure is used as fallback (pressure effect largely cancels).

class ModelHeightAdjustedSpecificHumidity : public ObsFunctionBase<float> {
 public:
    explicit ModelHeightAdjustedSpecificHumidity(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ModelHeightAdjustedSpecificHumidityParameters parameters_;
    ufo::Variables invars_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDSPECIFICHUMIDITY_H_
