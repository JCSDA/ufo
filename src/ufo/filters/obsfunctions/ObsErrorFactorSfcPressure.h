/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSFCPRESSURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSFCPRESSURE_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorSfcPressure ObsFunction
class ObsErrorFactorSfcPressureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorSfcPressureParameters, Parameters)

 public:
  /// the existence of min,max error values are required
  oops::RequiredParameter<float> error_min{"error_min", this};
  oops::RequiredParameter<float> error_max{"error_max", this};
  oops::Parameter<std::string> original_obserr{"original_obserr", "ObsErrorData", this};
  oops::Parameter<std::string> geovar_sfc_geomz{"geovar_sfc_geomz", "surface_altitude", this};
  oops::Parameter<std::string> geovar_temp{"geovar_temp", "virtual_temperature", this};
};

// -----------------------------------------------------------------------------

/// \brief Inflate the observation error for surface pressure as done by GSI-Observer.
///
/// This routine was designed to mimic the GSI observer code (i.e., setupps.f90) to inflate
/// the observation error for surface pressure using the following inputs:
///   Observed surface pressure, station height, and temperature (possibly missing).
///   Model first-guess fields interpolated to the observation location.
/// The starting obserror is then altered by this code with the "inflate error" action,
/// constrained by the values given for error_min and error_max (Pa).  For testing
/// purposes, the optional parameter of original_obserr group name such as ObsError to
/// override the default ObsErrorData can be used for tolerance check of reference results.
/// Internally, the code expects to find surface_altitude@GeoVaLs, however, some datasets
/// may have surface_geopotential_height@GeoVaLs in its place.  Also, the model vertical
/// profile of temperature is assumed to be virtual_temperature@GeoVaLs, but can be overriden
/// using geovar_temp option set to air_temperature for example.
///
/// ~~~~
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - filter: BlackList
///       filter variables:
///       - name: surface_pressure
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsErrorFactorSfcPressure@ObsFunction
///           options:
///             error_min: 100         # 1 mb
///             error_max: 300         # 3 mb
///             geovar_sfc_geomz: surface_geopotential_height     # default is surface_altitude
///
class ObsErrorFactorSfcPressure : public ObsFunctionBase {
 public:
  static const std::string classname() {return "ObsErrorFactorSfcPressure";}

  explicit ObsErrorFactorSfcPressure(const eckit::Configuration &config);
  ~ObsErrorFactorSfcPressure();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::unique_ptr<ObsErrorFactorSfcPressureParameters> options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSFCPRESSURE_H_
