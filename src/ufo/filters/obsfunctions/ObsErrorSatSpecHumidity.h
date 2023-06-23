/*
 * (C) Copyright 2023 NASA GMAO
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORSATSPECHUMIDITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORSATSPECHUMIDITY_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsErrorSatSpecHumidity ObsFunction
class ObsErrorSatSpecHumidityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorSatSpecHumidityParameters, Parameters)

 public:
  /// Inflate variables
  oops::RequiredParameter<std::string> inflatevars{"variable", this};
  oops::RequiredParameter<std::string> inputerr_name{"input_error_name", this};
  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};
};

// -----------------------------------------------------------------------------

/// \brief Assign observation error for specific humidity as done by GSI-Observer.
///
/// ObsErrorSatSpecHumidity assigns the errors with an inflation using
/// saturated specific humdity. This field is obtained through GSI geovals as a
/// profile
///
/// Notes:
/// - This obs function is designed for specific humidity observations only. When
///   used in a filter, please make sure "filter variables" only contains the
///   variable name for specific humidity.
///
/// ### example configurations for using this obs function in a filter: ###
///   - filter: Perform Action
///     filter variables:
///     - name: specificHumidity
///     action:
///       name: assign error
///       error function:
///         name: ObsFunction/ObsErrorSaturatedSpecificHumidity
///         options:
///           variable: specificHumidity
///           input_error_name: GsiInputObsError
///
class ObsErrorSatSpecHumidity : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorSatSpecHumidity";}

  explicit ObsErrorSatSpecHumidity(const eckit::Configuration &config);
  ~ObsErrorSatSpecHumidity();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::unique_ptr<ObsErrorSatSpecHumidityParameters> options_;
  bool isAscending_ = true;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORSATSPECHUMIDITY_H_
