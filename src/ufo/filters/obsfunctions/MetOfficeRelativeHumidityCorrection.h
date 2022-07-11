/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_METOFFICERELATIVEHUMIDITYCORRECTION_H_
#define UFO_FILTERS_OBSFUNCTIONS_METOFFICERELATIVEHUMIDITYCORRECTION_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the MetOfficeRelativeHumidityCorrection ObsFunction
class MetOfficeRelativeHumidityCorrectionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MetOfficeRelativeHumidityCorrectionParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> model_pressure
    {"model pressure",
     "Name of model pressure.",
     this};

  oops::RequiredParameter<std::string> model_specific_humidity
    {"model specific humidity",
     "Name of model specific humidity.",
     this};

  oops::RequiredParameter<std::string> model_relative_humidity
    {"model relative humidity",
     "Name of model relative humidity.",
     this};

  oops::RequiredParameter<std::string> model_temperature
    {"model temperature",
     "Name of model temperature.",
     this};

  oops::RequiredParameter<std::string> observed_pressure
    {"observed pressure",
     "Name of observed pressure.",
     this};

  oops::Parameter<bool> capsupersat
    {"capsupersat",
     "Cap relative humidity at 100%",
     false,
     this};
};

// -----------------------------------------------------------------------------

/// \brief Compute a correction factor that accounts for differences in relative humidity H(x)
/// produced by the UM interface and the Met Office OPS system.
///
/// There are differences in values of relative humidity H(x) produced by the UM interface and
/// the Met Office OPS system. The differences are caused by the order in which
/// (1) RH is computed from specific humidity (2) temporal and spatial interpolation are
/// performed. The computation of relative humidity from specific humidity is nonlinear,
/// which can lead to differences in H(x) of up to 5%.
///
/// This ObsFunction computes two values of relative humidity H(x). The first reproduces
/// what occurs in the UM interface, i.e. it vertically interpolates the relative humidity
/// GeoVaL at each location. The second reproduces what occurs in OPS, i.e. it computes
/// relative humidity from GeoVaLs of specific humidity, temperature, and pressure and
/// then performs vertical interpolation. The output of the ObsFunction is the difference
/// between the two interpolated H(x) values.
///
/// The H(x) difference can be added to observed relative humidity (using a LinearCombination
/// ObsFunction) prior to running QC filters that rely on relative humidity O-B such as the
/// Background Check. After those filters have run, the H(x) difference can be subtracted back off.
class MetOfficeRelativeHumidityCorrection : public ObsFunctionBase<float> {
 public:
  explicit MetOfficeRelativeHumidityCorrection(const eckit::LocalConfiguration &);
  ~MetOfficeRelativeHumidityCorrection();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  MetOfficeRelativeHumidityCorrectionParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_METOFFICERELATIVEHUMIDITYCORRECTION_H_
