/*
 * (C) Copyright 2022 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELHUMIDITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELHUMIDITY_H_

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

/// \brief Options controlling ObsErrorModelHumidity ObsFunction
class ObsErrorModelHumidityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorModelHumidityParameters, Parameters)

 public:
  /// the name of the xvar
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  /// vector of X-steps
  oops::RequiredParameter<std::vector<float>> xvals{"xvals", this};
  /// vector of error values corresponding to vector of X-steps
  oops::RequiredParameter<std::vector<float>> errors{"errors", this};
  /// Method used for saturated vapor pressure [Optional]
  oops::Parameter<std::string> Method{"Method", "default", this};
  /// Formulation used for saturated vapor pressure [Optional]
  oops::Parameter<std::string> Formulation{"Formulation", "", this};
};

// -----------------------------------------------------------------------------

/// \brief Assign observation error for specific humidity as done by GSI-Observer.
///
/// This routine was designed to mimic the way GSI observer code (i.e., setupq.f90)
/// assigns observation error for specific humidity observations in two steps.
/// The first step is to interpolate errors from GSI fix-file prepobs_errtable.txt
/// to observation pressure with parameters and routines borrowed from
/// ObsFunction/ObsErrorModelStepwiseLinear. The second step is to scale the errors
/// with saturation specific humidity estimated from model fields because humidity
/// errors from prepobs_errtable.txt are for relative humidity observations.
/// Forecast temperature and pressure fields are required to estimate saturation
/// specific humidity.
///
/// Notes:
/// - This obs function is designed for specific humidity observations only. When
///   used in a filter, please make sure "filter variables" only contains the
///   variable name for specific humidity.
/// - Required parameters (xvar, xvals, and errors) and routines utilizing them are
///   borrowed from ObsFunction/ObsErrorModelStepwiseLinear, so YAML configuration
///   for this function must also follow that of ObsErrorModelStepwiseLinear. One
///   modification is to allow size-1 lists of xvals and errors to accommodate the
///   situation of fixed input error.
/// - The method or formulation to compute saturation specific humidity is
///   configurable via YAML. If not specified, the default formula will be used.
///
/// ### example configurations for using this obs function in a filter: ###
///   - filter: Perform Action
///     filter variables:
///     - name: specificHumidity
///     action:
///       name: assign error
///       error function:
///         name: ObsFunction/ObsErrorModelHumidity
///         options:
///           xvar:
///             name: MetaData/pressure
///           xvals: [85000, 50000, 25000]   #Pressure (Pa)
///           errors: [0.18, 0.19, 0.2]      #RH error
///           Method: UKMO

class ObsErrorModelHumidity : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorModelHumidity";}

  explicit ObsErrorModelHumidity(const eckit::Configuration &config);
  ~ObsErrorModelHumidity();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsErrorModelHumidityParameters options_;
  bool isAscending_ = true;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORMODELHUMIDITY_H_
