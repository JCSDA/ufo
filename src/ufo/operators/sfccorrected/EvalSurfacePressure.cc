/*
 * (C) Copyright 2025, UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "ufo/operators/sfccorrected/EvalSurfacePressure.h"

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

namespace {
SurfaceOperatorMaker<stationPressure_WRFDA> makerSP_WRFDA_("stationPressure_WRFDA");
SurfaceOperatorMaker<stationPressure_UKMO> makerSP_UKMO_("stationPressure_UKMO");
SurfaceOperatorMaker<stationPressure_GSL> makerSP_GSL_("stationPressure_GSL");

// Returns true when the geovariable name indicates a geopotential height,
// in which case heights need to be converted to geometric height before use.
bool isGeopotentialHeightVariable(const std::string & geovar) {
  return geovar.find("geopotential") != std::string::npos;
}

// Computes virtual temperature from air temperature and specific humidity at each location.
// Locations with a missing air temperature or specific humidity are left unchanged.
void computeVirtualTemperature(const std::vector<float> & airTemp,
                               const std::vector<float> & specificHumidity,
                               std::vector<float> & virtualTemp) {
  const float missing = util::missingValue<float>();
  for (size_t iloc = 0; iloc < virtualTemp.size(); ++iloc) {
    if (airTemp[iloc] != missing && specificHumidity[iloc] != missing) {
      virtualTemp[iloc] = formulas::VirtualTemp_From_Sh_AT(specificHumidity[iloc], airTemp[iloc]);
    }
  }
}
}  // namespace

// ----------------------------------------
// Pressure operator using WRFDA method
// ----------------------------------------
stationPressure_WRFDA::stationPressure_WRFDA(const std::string & name,
                                             const Parameters_ & params)
     : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_pressure"));
  vars.push_back(oops::Variable("virtual_temperature"));
  vars.push_back(oops::Variable("air_pressure_at_surface"));
  requiredVars_ += vars;
}

// --------------------------------------------------------------------------------
// Before calculating Hofx, Conduct terrain height correction for surface pressure
//
// Date: June 2019: Created
// Method: hydrosatic equation
//
//  adjusted_model_surface_height_pressure = obs_pressure * exp [-grav/rd *
//   (model_height_surface - model_height_level1) / (model_virtual_T + obs_virtual_T)/2]
//
//  Where:
//  model_height_surface = model surface height
//  model_height_level1 = station height
//  model_virtual_T = virtual temperature at model surface height
//  obs_virtual_T = virtual temperature at station height
//  adjusted_model_surface_height_pressure =
//       pressure interpolated from station height to model surface height
//  obs_pressure = pressure at station height
//  grav = gravitational acceleration
//  rd = gas constant per mole
//
// --------------------------------------------------------------------------------
void stationPressure_WRFDA::simobs(const ufo::GeoVaLs & gv,
                                   const ioda::ObsSpace & obsdb,
                                   std::vector<float> & hofx) const {
  oops::Log::trace() << "stationPressure_WRFDA::simobs start" << std::endl;

  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const float missingFloat = util::missingValue<float>();

  // Create arrays needed
  std::vector<float> model_pressure_surface(nobs), model_height_level1(nobs),
                     model_height_surface(nobs), model_virtual_T(nobs);
  std::vector<float> obs_pressure(nobs), obs_virtual_T(nobs),
                     obs_lats(nobs), obs_height(nobs);
  std::vector<float> adjusted_model_surface_virtual_T(nobs);
  std::vector<float> adjusted_station_height_virtual_T(nobs);
  std::vector<float> adjusted_model_surface_height_pressure(nobs);

  this->getDataValues(gv, obsdb, params_, obs_lats, obs_height, obs_virtual_T,
                      obs_pressure, model_height_level1, model_height_surface,
                      model_pressure_surface, model_virtual_T);
  // Loop to calculate hofx
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = missingFloat;
    if (obs_height[iloc] != missingFloat && model_pressure_surface[iloc] != missingFloat &&
         model_height_surface[iloc] != missingFloat &&
         model_height_level1[iloc] != missingFloat) {
      // Find model surface virtual temperature
      adjusted_model_surface_virtual_T[iloc] = model_virtual_T[iloc] +
               ufo::Constants::Lclr * (model_height_level1[iloc] - model_height_surface[iloc]);
      if (obs_virtual_T[iloc] != missingFloat) {
        adjusted_station_height_virtual_T[iloc] = 0.5 * (
               adjusted_model_surface_virtual_T[iloc] + obs_virtual_T[iloc]);
      } else {
        adjusted_station_height_virtual_T[iloc] = adjusted_model_surface_virtual_T[iloc];
      }
      // Extrapolate pressure from station height to model surface height
      if (obs_pressure[iloc] != missingFloat) {
        float val = (model_height_surface[iloc] - obs_height[iloc]) *
                (ufo::Constants::grav /
                (ufo::Constants::rd * adjusted_station_height_virtual_T[iloc]));
        adjusted_model_surface_height_pressure[iloc] = obs_pressure[iloc] * std::exp(-val);
      } else {
        adjusted_model_surface_height_pressure[iloc] = missingFloat;
      }
      if (adjusted_model_surface_height_pressure[iloc] != missingFloat) {
        hofx[iloc] = obs_pressure[iloc] - adjusted_model_surface_height_pressure[iloc] +
                     model_pressure_surface[iloc];
      } else {
        hofx[iloc] = model_pressure_surface[iloc];
      }
    }
  }
  oops::Log::trace() << "stationPressure_WRFDA::simobs done" << std::endl;
}

// --------------------------------------------------------------------------------
void stationPressure_WRFDA::getTemperatureValues(const size_t nobs,
                                                  const int surface_level_index,
                                                  const ufo::GeoVaLs & gv,
                                                  const ioda::ObsSpace & obsdb,
                                                  std::vector<float> & modelVirtualTemp,
                                                  std::vector<float> & obsVirtualTemp) const {
  gv.getAtLevel(modelVirtualTemp, oops::Variable("virtual_temperature"), surface_level_index);

  // based on what variables ObsSpace contains, obsVirtualTemp value is assigned
  enum class TempSource { ObsVirtualT, ObsAirTAndSpecificHumidity, ObsAirT, ModelVirtualT };
  const TempSource source =
      obsdb.has("ObsValue", "virtualTemperatureAt2M")
          ? TempSource::ObsVirtualT
      : !obsdb.has("ObsValue", "airTemperatureAt2M")
          ? TempSource::ModelVirtualT
      : obsdb.has("ObsValue", "specificHumidityAt2M")
          ? TempSource::ObsAirTAndSpecificHumidity
          : TempSource::ObsAirT;

  const float missing = util::missingValue<float>();
  obsVirtualTemp.assign(nobs, missing);

  switch (source) {
    case TempSource::ObsVirtualT: {
      obsdb.get_db("ObsValue", "virtualTemperatureAt2M", obsVirtualTemp);
      break;
    }
    case TempSource::ObsAirTAndSpecificHumidity: {
      // Calculate virtual temperature from air temperature and specific humidity
      std::vector<float> obsAirTemp(nobs, missing);
      std::vector<float> obsHumidity(nobs, missing);
      obsdb.get_db("ObsValue", "airTemperatureAt2M", obsAirTemp);
      obsdb.get_db("ObsValue", "specificHumidityAt2M", obsHumidity);
      computeVirtualTemperature(obsAirTemp, obsHumidity, obsVirtualTemp);
      break;
    }
    case TempSource::ObsAirT: {
      // Specific humidity not available; use air temperature as virtual temperature.
      obsdb.get_db("ObsValue", "airTemperatureAt2M", obsVirtualTemp);
      break;
    }
    case TempSource::ModelVirtualT: {
      oops::Log::warning() << "stationPressure_WRFDA::getDataValues virtualTemperature and "
                           << "airTemperature not found, set to model virtualTemperature"
                           << std::endl;
      obsVirtualTemp = modelVirtualTemp;
      break;
    }
  }
}

// --------------------------------------------------------------------------------
void stationPressure_WRFDA::getDataValues(const ufo::GeoVaLs & gv,
                                          const ioda::ObsSpace & obsdb,
                                          const Parameters_ & params_,
                                          std::vector<float> & obsLats,
                                          std::vector<float> & obsHeight,
                                          std::vector<float> & obsVirtualTemp,
                                          std::vector<float> & obsPressure,
                                          std::vector<float> & modelHeightLevel1,
                                          std::vector<float> & modelHeightSurface,
                                          std::vector<float> & modelPressureSurface,
                                          std::vector<float> & modelVirtualTemp) const {
  const size_t nobs = obsdb.nlocs();
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  const int surface_level_index = gv.nlevs(geomz_var) - 1;

  obsdb.get_db("MetaData", "stationElevation", obsHeight);
  obsdb.get_db("MetaData", "latitude", obsLats);
  obsdb.get_db("ObsValue", "stationPressure", obsPressure);

  // Get level 1 height.  If geopotential then convert to geometric height.
  gv.getAtLevel(modelHeightLevel1, geomz_var, surface_level_index);
  if (isGeopotentialHeightVariable(params_.geovarGeomZ.value())) {
    oops::Log::trace()  <<
        "ObsSfcCorrected::simulateObs_WRFDA do geopotential conversion for model level 1"
        << std::endl;
    convertToGeometricHeight(obsLats, modelHeightLevel1);
  }
  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(modelHeightSurface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs_WRFDA do geopotential conversion for "
                        << "model surface level" << std::endl;
    convertToGeometricHeight(obsLats, modelHeightSurface);
  }
  gv.get(modelPressureSurface, oops::Variable("air_pressure_at_surface"));
  this->getTemperatureValues(nobs, surface_level_index, gv, obsdb,
                             modelVirtualTemp, obsVirtualTemp);
}

// ----------------------------------------
// Pressure operator using UKMO method
// ----------------------------------------
stationPressure_UKMO::stationPressure_UKMO(const std::string & name,
                                           const Parameters_ & params)
     : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_pressure"));
  vars.push_back(oops::Variable("virtual_temperature"));
  vars.push_back(oops::Variable("air_pressure_at_surface"));
  requiredVars_ += vars;
}

// --------------------------------------------------------------------------------
// Before calculating Hofx, Conduct terrain height correction for surface pressure
//
// Reference: Ingleby,2013. UKMO Technical Report No: 582. Appendix 1.
//
// Method: integrate the hydrosatic equation dp/dz=-rho*g/RT
// to get P_m2o first, equation:
//
//  (model_P_obs_height/model_pressure_surface) =
//         (model_T_obs_height/model_T_surface)** (grav/rd*L)
//
//  Where:
//  model_P_obs_height = model surface pressure at station height
//  model_pressure_surface = model surface pressure
//  model_T_surface = temperature at model surface height; derived at 2000m height
//  model_T_obs_height = model surface temperature at station height
//  grav  = gravitational acceleration
//  rd    = gas constant per mole
//  Lclr  = constant lapse rate (0.0065 K/m)
//
//  To avoid dirunal/local variations, use TV_2000 (2000 m above the model surface height)
//  instead of direct model_T_surface
//
//  model_T_surface = adjusted_model_virtual_T_2000m *
//      (model_pressure_surface / adjusted_model_pressure_2000m) ** (rd*L/grav)
//
// Where:
//  adjusted_model_pressure_2000m = background pressure at 2000 m
//  adjusted_model_virtual_T_2000m = background virtual temperature at 2000 m
//  obs_pressure =  pressure at station height
//
//  Finally, in practice, adjust obs_pressure to the model surface height using
//
//  adjusted_model_surface_height_pressure =
//          obs_pressure * (model_pressure_surface / model_P_obs_height)
//
// --------------------------------------------------------------------------------
void stationPressure_UKMO::simobs(const ufo::GeoVaLs & gv,
                                  const ioda::ObsSpace & obsdb,
                                  std::vector<float> & hofx) const {
  oops::Log::trace() << "stationPressure_UKMO::simobs start" << std::endl;

  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  const size_t model_nlevs = gv.nlevs(geomz_var);
  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  const oops::Variable model_p_var = oops::Variable("air_pressure");
  const oops::Variable model_virtual_T_var = oops::Variable("virtual_temperature");
  const oops::Variable model_height_var = oops::Variable(params_.geovarGeomZ.value());

  // Create arrays needed
  std::vector<float> obs_lats(nobs), obs_height(nobs), obs_pressure(nobs);
  std::vector<float> adjusted_model_surface_height_pressure(nobs, missingFloat);
  std::vector<float> model_height_surface(nobs, missingFloat);
  std::vector<float> model_pressure_surface(nobs, missingFloat);
  std::vector<double> adjusted_model_pressure_2000m(nobs, missingDouble);
  std::vector<double> adjusted_model_virtual_T_2000m(nobs, missingDouble);

  this->getDataValues(gv, obsdb, params_, obs_lats, obs_height, obs_pressure,
                      model_height_surface, model_pressure_surface);
  // Decide if geopotential -> geometric conversion needed
  const bool convertLevel1GeopotentialHeight =
      isGeopotentialHeightVariable(params_.geovarGeomZ.value());

  const double height_used = 2000.0;
  int index = 0;
  double weight = 0.0;
  double power_exponent = ( (ufo::Constants::Lclr * ufo::Constants::rd) / ufo::Constants::grav);
  std::vector<double> model_pressure(model_nlevs), model_virtual_T(model_nlevs),
                      model_height_level1(model_nlevs);

  // Loop to calculate hofx
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = model_pressure_surface[iloc];
    if (obs_height[iloc] != missingFloat && model_pressure_surface[iloc] != missingFloat &&
        model_height_surface[iloc] != missingFloat && obs_pressure[iloc] != missingFloat) {
      // Get model data at this location
      gv.getAtLocation(model_height_level1, model_height_var, iloc);
      gv.getAtLocation(model_pressure, model_p_var, iloc);
      gv.getAtLocation(model_virtual_T, model_virtual_T_var, iloc);

      // Convert geopotential to geometric height if needed
      if (convertLevel1GeopotentialHeight) {
        convertToGeometricHeight(obs_lats[iloc], model_height_level1);
      }
      vert_interp_weights_f90(model_nlevs, height_used, model_height_level1.data(),
                              index, weight);
      // Vertical interpolation to get model pressure and temperature at 2000 m
      vert_interp_apply_f90(model_nlevs, model_pressure.data(),
                            adjusted_model_pressure_2000m[iloc], index, weight);
      vert_interp_apply_f90(model_nlevs, model_virtual_T.data(),
                            adjusted_model_virtual_T_2000m[iloc], index, weight);

      // Calculate background temperature at model surface height and at station height
      if (adjusted_model_pressure_2000m[iloc] != missingDouble &&
          adjusted_model_virtual_T_2000m[iloc] != missingDouble) {
        // Model Temperature at Surface, and Observed height
        double model_T_surface = adjusted_model_virtual_T_2000m[iloc] *
          std::pow((model_pressure_surface[iloc] / adjusted_model_pressure_2000m[iloc]),
                    power_exponent);
        double model_T_obs_height = model_T_surface + (ufo::Constants::Lclr *
                          (model_height_surface[iloc] - obs_height[iloc]));

        // Model Pressure at Observed height
        double power_exponent_inverse = 1.0 / power_exponent;
        double model_P_obs_height = model_pressure_surface[iloc] *
                std::pow((model_T_obs_height / model_T_surface), power_exponent_inverse);
        adjusted_model_surface_height_pressure[iloc] = (obs_pressure[iloc] *
                model_pressure_surface[iloc]) / model_P_obs_height;
      }
    }
    if (adjusted_model_surface_height_pressure[iloc] != missingFloat) {
      hofx[iloc] = obs_pressure[iloc] - adjusted_model_surface_height_pressure[iloc] +
                   model_pressure_surface[iloc];
    }
  }
  oops::Log::trace() << "stationPressure_UKMO::simobs done" << std::endl;
}

// --------------------------------------------------------------------------------
void stationPressure_UKMO::getDataValues(const ufo::GeoVaLs & gv,
                                          const ioda::ObsSpace & obsdb,
                                          const Parameters_ & params_,
                                          std::vector<float> & obsLats,
                                          std::vector<float> & obsHeight,
                                          std::vector<float> & obsPressure,
                                          std::vector<float> & modelHeightSurface,
                                          std::vector<float> & modelPressureSurface) const {
  obsdb.get_db("MetaData", "stationElevation", obsHeight);
  obsdb.get_db("MetaData", "latitude", obsLats);
  obsdb.get_db("ObsValue", "stationPressure", obsPressure);

  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(modelHeightSurface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs_UKMO do geopotential conversion for "
                        << "model surface level" << std::endl;
    convertToGeometricHeight(obsLats, modelHeightSurface);
  }
  gv.get(modelPressureSurface, oops::Variable("air_pressure_at_surface"));
}

// ----------------------------------------
// Pressure operator using GSL method
// ----------------------------------------
stationPressure_GSL::stationPressure_GSL(const std::string & name,
                                         const Parameters_ & params)
     : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_pressure"));
  vars.push_back(oops::Variable("virtual_temperature"));
  vars.push_back(oops::Variable("air_pressure_at_surface"));
  requiredVars_ += vars;
}

// --------------------------------------------------------------------------------
// Before calculating Hofx, Conduct terrain height correction for surface pressure
//
//  adjusted_model_surface_height_pressure = exp(log(model_pressure_surface) -
//       ((obs_height - model_height_surface) * (grav/rd) / avg_virtual_T))
//
//  Where:
//  model_height_surface = model surface height
//  obs_height = observed station height
//  model_pressure_surface = model surface pressure
//  avg_virtual_T = average of modeled/observed virtual temperature
//  grav = gravitational acceleration
//  rd = gas constant per mole
//  adjusted_model_surface_height_pressure = model surface pressure at station height
//
// --------------------------------------------------------------------------------
void stationPressure_GSL::simobs(const ufo::GeoVaLs & gv,
                                 const ioda::ObsSpace & obsdb,
                                 std::vector<float> & hofx) const {
  oops::Log::trace() << "stationPressure_GSL::simobs start" << std::endl;

  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  const int model_nlevs = gv.nlevs(geomz_var);
  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  // Create arrays needed
  std::vector<double> model_pressure_surface(nobs), model_virtual_T(model_nlevs);
  std::vector<float> obs_pressure(nobs), obs_T(nobs), obs_virtual_T(nobs),
                     obs_lats(nobs), obs_height(nobs), model_virtual_T_surface(nobs),
                     model_height_surface(nobs);

  this->getDataValues(gv, obsdb, params_, obs_lats, obs_height, obs_pressure, obs_virtual_T,
              obs_T, model_height_surface, model_pressure_surface, model_virtual_T_surface);

  int index = 0;
  double weight = 0.0;
  std::vector<double> model_T_interpolated(nobs);
  std::vector<double> model_log_pressure(model_nlevs), model_pressure(model_nlevs);
  std::vector<float> adjusted_model_surface_height_pressure(nobs), avg_virtual_T(nobs);

  // Loop to calculate corrections
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    gv.getAtLocation(model_pressure, oops::Variable("air_pressure"), iloc);
    const bool model_pressure_missing = std::any_of(model_pressure.begin(), model_pressure.end(),
                      [missingDouble](double p) { return p == missingDouble; });
    if (obs_virtual_T[iloc] == missingFloat && obs_T[iloc] != missingFloat) {
      obs_virtual_T[iloc] = obs_T[iloc];
    }
    gv.getAtLocation(model_virtual_T, oops::Variable("virtual_temperature"), iloc);
    avg_virtual_T[iloc] = model_virtual_T_surface[iloc];

    if (obs_virtual_T[iloc] != missingFloat && obs_virtual_T[iloc] > 150.0
       && obs_virtual_T[iloc] < 350.0) {
      avg_virtual_T[iloc] = (model_virtual_T_surface[iloc] + obs_virtual_T[iloc]) / 2.0;
    } else {
      if (obs_pressure[iloc] != missingFloat && obs_height[iloc] != missingFloat &&
          model_height_surface[iloc] != missingFloat && !model_pressure_missing) {
        std::transform(model_pressure.begin(), model_pressure.end(),
            model_log_pressure.begin(), [](double p) { return std::log(p); });
        vert_interp_weights_f90(model_nlevs, std::log(static_cast<double>(obs_pressure[iloc])),
                 model_log_pressure.data(), index, weight);
        // Extrapolate model temperature to obs_height
        vert_interp_apply_f90(model_nlevs, model_virtual_T.data(),
                              model_T_interpolated[iloc], index, weight);
        avg_virtual_T[iloc] = (model_virtual_T_surface[iloc] + model_T_interpolated[iloc]) / 2.0;

        if (obs_height[iloc] < model_height_surface[iloc]) {
          // Extrapolate to surface if observation is below lowest model layer
          avg_virtual_T[iloc] = avg_virtual_T[iloc] - ((ufo::Constants::Lclr / 2.0) *
               (obs_height[iloc] - model_height_surface[iloc]));
        }
      }
    }
  }
  // Loop to calculate hofx
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = missingFloat;
    if (obs_height[iloc] != missingFloat && model_pressure_surface[iloc] != missingDouble &&
         model_height_surface[iloc] != missingFloat && avg_virtual_T[iloc] != 0.0) {
      // Extrapolate pressure from model surface to observation height
      adjusted_model_surface_height_pressure[iloc] = std::exp(
         std::log(model_pressure_surface[iloc]) -
         ((ufo::Constants::grav / ufo::Constants::rd) *
          (obs_height[iloc] - model_height_surface[iloc])) /
           avg_virtual_T[iloc]);
    } else {
      adjusted_model_surface_height_pressure[iloc] = model_pressure_surface[iloc];
    }
    hofx[iloc] = adjusted_model_surface_height_pressure[iloc];
  }
  oops::Log::trace() << "stationPressure_GSL::simobs done" << std::endl;
}

// --------------------------------------------------------------------------------
void stationPressure_GSL::getDataValues(const ufo::GeoVaLs & gv,
                                        const ioda::ObsSpace & obsdb,
                                        const Parameters_ & params_,
                                        std::vector<float> & obsLats,
                                        std::vector<float> & obsHeight,
                                        std::vector<float> & obsPressure,
                                        std::vector<float> & obsVirtualTemp,
                                        std::vector<float> & obsTemp,
                                        std::vector<float> & modelHeightSurface,
                                        std::vector<double> & modelPressureSurface,
                                        std::vector<float> & modelVirtualTempSurface) const {
  const size_t nobs = obsdb.nlocs();
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  const int surface_level_index = gv.nlevs(geomz_var) - 1;
  const float missing = util::missingValue<float>();

  obsdb.get_db("MetaData", params_.obsHeightName.value(), obsHeight);
  obsdb.get_db("MetaData", "latitude", obsLats);
  obsdb.get_db("ObsValue", "stationPressure", obsPressure);

  obsVirtualTemp.assign(nobs, missing);
  obsTemp.assign(nobs, missing);
  const ufo::Variable & tempVar = params_.obsTemperature.value();
  const ufo::Variable & humidityVar = params_.obsSpecificHumidity.value();
  if (obsdb.has(tempVar.group(), tempVar.variable())) {
    obsdb.get_db(tempVar.group(), tempVar.variable(), obsTemp);
  }
  // Derive observed virtual temperature if specific humidity is available.
  if (obsdb.has(humidityVar.group(), humidityVar.variable())) {
    std::vector<float> obsHumidity(nobs, missing);
    obsdb.get_db(humidityVar.group(), humidityVar.variable(), obsHumidity);
    computeVirtualTemperature(obsTemp, obsHumidity, obsVirtualTemp);
  }

  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(modelHeightSurface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs_GSL do geopotential conversion for "
                        << "model surface level" << std::endl;
    convertToGeometricHeight(obsLats, modelHeightSurface);
  }
  gv.get(modelPressureSurface, oops::Variable("air_pressure_at_surface"));
  gv.getAtLevel(modelVirtualTempSurface, oops::Variable("virtual_temperature"),
                surface_level_index);
}

}  // namespace ufo
