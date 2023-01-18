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

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorSfcPressure ObsFunction
class ObsErrorFactorSfcPressureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorSfcPressureParameters, Parameters)

 public:
  oops::Parameter<std::string> original_obserr{"original_obserr", "ObsErrorData", this};
  oops::Parameter<std::string> geovar_geomz{"geovar_geomz", "height", this};
  oops::Parameter<std::string> geovar_sfc_geomz{"geovar_sfc_geomz", "surface_altitude", this};
};

// -----------------------------------------------------------------------------

/// \brief Inflate the observation error for surface pressure as done by GSI-Observer.
///
/// This routine was designed to mimic the GSI observer code (i.e., setupps.f90) to inflate
/// the observation error for surface pressure using the following inputs:
///   Observed surface pressure, station height, virtual temperature or air temperature.
///   Model first-guess fields interpolated to the observation location.
/// The starting obserror is then altered by this code with the "inflate error" action,
/// For testing purposes, the optional parameter of original_obserr group name such as ObsError
/// to override the default ObsErrorData can be used for tolerance check of reference results.
/// Internally, the code expects to find GeoVals/surface_altitude and GeoVaLs/height, however,
/// some datasets may have GeoVaLs/surface_geopotential_height and GeoVaLs/geopotential_height
/// in its place.
///
/// ~~~~
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - filter: BlackList
///       filter variables:
///       - name: stationPressure
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsFunction/ObsErrorFactorSfcPressure
///           options:
///             geovar_geomz: geopotential_height     # default is height
///             geovar_sfc_geomz: surface_geopotential_height     # default is surface_altitude
///
class ObsErrorFactorSfcPressure : public ObsFunctionBase<float> {
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
