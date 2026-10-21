/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "ufo/operators/sfccorrected/EvalSurfaceTemperature.h"

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

namespace {
SurfaceOperatorMaker<airTemperatureAt2M_WRFDA> makerT2M_WRFDA_("airTemperatureAt2M_WRFDA");
SurfaceOperatorMaker<airTemperatureAt2M_UKMO> makerT2M_UKMO_("airTemperatureAt2M_UKMO");
SurfaceOperatorMaker<airTemperatureAt2M_GSL> makerT2M_GSL_("airTemperatureAt2M_GSL");

// Returns true when the geovariable name indicates a geopotential height,
// in which case heights need to be converted to geometric height before use.
bool isGeopotentialHeightVariable(const std::string & geovar) {
  return geovar.find("geopotential") != std::string::npos;
}

/// Result of the GSL "Local" lapse-rate calculation, shared by simobs and settraj.
struct GslLocalLapseRate {
  int lapseTopLevelIndex;
  std::vector<float> lapseRate;   // (T_lev1 - T_top) / (z_top - z_lev1) in K/m, clamped
                                  // if requested; missing where inputs are invalid
  std::vector<float> lapseLayerDepth;  // z_top - z_lev1 in m (meaningful where lapseRate is valid)
  std::vector<bool> lapseRateClamped;  // true where the threshold clamp changed the lapse rate
};

/// Lapse rate between the lowest model level and the level
/// "temperature local lapse rate level" layers above it.
GslLocalLapseRate computeGslLocalLapseRate(const GeoVaLs & gv,
                                           const std::vector<float> & lats,
                                           const std::string & geovarGeomZ,
                                           const GslCorrectionParameters & gsl,
                                           const int surfaceLevelIndex,
                                           const std::string & caller) {
  const size_t nobs = lats.size();
  const float missing = util::missingValue<float>();

  const int local_highest_level = gsl.temperatureLocalLapseRateLevel.value();
  if (surfaceLevelIndex <= local_highest_level || local_highest_level < 1) {
    throw eckit::AssertionFailed(caller + " local_highest_level must be in range of 1,"
                                 " surface_level_index -1", Here());
  }

  GslLocalLapseRate result;
  result.lapseTopLevelIndex = surfaceLevelIndex - local_highest_level + 1;
  result.lapseRate.assign(nobs, missing);
  result.lapseLayerDepth.assign(nobs, 0.0f);
  result.lapseRateClamped.assign(nobs, false);

  // Heights of the two levels (geopotential -> geometric if needed).
  const oops::Variable geomz_var = oops::Variable(geovarGeomZ);
  std::vector<float> model_height_level1(nobs);
  std::vector<float> model_height_toplayer(nobs);
  gv.getAtLevel(model_height_level1, geomz_var, surfaceLevelIndex);
  gv.getAtLevel(model_height_toplayer, geomz_var, result.lapseTopLevelIndex);
  if (isGeopotentialHeightVariable(geovarGeomZ)) {
    convertToGeometricHeight(lats, model_height_level1);
    convertToGeometricHeight(lats, model_height_toplayer);
  }

  // Temperatures of the two levels.
  const oops::Variable model_T_var = oops::Variable("air_temperature");
  std::vector<float> model_T_level1(nobs);
  std::vector<float> model_T_toplayer(nobs);
  gv.getAtLevel(model_T_level1, model_T_var, surfaceLevelIndex);
  gv.getAtLevel(model_T_toplayer, model_T_var, result.lapseTopLevelIndex);

  const bool use_threshold = gsl.temperatureLapseRateThreshold.value();
  // K/km -> K/m
  const float minthresh = static_cast<float>(gsl.minThreshold.value() / ufo::Constants::thousand);
  const float maxthresh = static_cast<float>(gsl.maxThreshold.value() / ufo::Constants::thousand);
  for (size_t iobs = 0; iobs < nobs; ++iobs) {
    if (model_height_toplayer[iobs] == missing || model_height_level1[iobs] == missing ||
        model_T_toplayer[iobs] == missing || model_T_level1[iobs] == missing ||
        model_height_toplayer[iobs] == model_height_level1[iobs]) {
      continue;
    }
    const float dz = model_height_toplayer[iobs] - model_height_level1[iobs];
    float lr = (model_T_level1[iobs] - model_T_toplayer[iobs]) / dz;
    if (use_threshold) {
      const float lr_clamped = std::clamp(lr, minthresh, maxthresh);
      result.lapseRateClamped[iobs] = (lr_clamped != lr);
      lr = lr_clamped;
    }
    result.lapseRate[iobs] = lr;
    result.lapseLayerDepth[iobs] = dz;
  }
  return result;
}
}  // namespace

// ----------------------------------------
// Temperature operator using WRFDA method
// ----------------------------------------
airTemperatureAt2M_WRFDA::airTemperatureAt2M_WRFDA(const std::string & name,
                                                   const Parameters_ & params)
    : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_temperature"));
  requiredVars_ += vars;
}

void airTemperatureAt2M_WRFDA::simobs(const ufo::GeoVaLs & gv,
                                      const ioda::ObsSpace & obsdb,
                                      std::vector<float> & hofx) const {
  oops::Log::trace() << "airTemperatureAt2M_WRFDA::simobs start" << std::endl;

  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const float missingFloat = util::missingValue<float>();

  std::vector<float> obs_lats(nobs);
  obsdb.get_db("MetaData", "latitude", obs_lats);

  // Create arrays needed
  std::vector<float> model_height_level1(nobs), model_height_surface(nobs),
                     model_T_level1(nobs), obs_height(nobs);

  // Get level 1 height.  If geopotential then convert to geometric height.
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  const int surface_level_index = gv.nlevs(geomz_var) - 1;
  gv.getAtLevel(model_height_level1, geomz_var, surface_level_index);
  if (isGeopotentialHeightVariable(params_.geovarGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs_WRFDA do geopotential"
                        << " conversion for model level 1" << std::endl;
    convertToGeometricHeight(obs_lats, model_height_level1);
  }

  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(model_height_surface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs_WRFDA do geopotential"
                        << " conversion for model surface level" << std::endl;
    convertToGeometricHeight(obs_lats, model_height_surface);
  }

  // Read other data in
  oops::Variable model_T_var = oops::Variable("air_temperature");
  gv.getAtLevel(model_T_level1, model_T_var, surface_level_index);
  obsdb.get_db("MetaData", params_.obsHeightName.value(), obs_height);

  // Loop to calculate hofx
  std::vector<float> model_T_surface(nobs);
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = missingFloat;
    if (obs_height[iloc] != missingFloat && model_T_level1[iloc] != missingFloat &&
        model_height_level1[iloc] != missingFloat && model_height_surface[iloc] != missingFloat) {
      // Find model surface T using lowest T in model temperature profile
      model_T_surface[iloc] = model_T_level1[iloc] +
              ufo::Constants::Lclr * (model_height_level1[iloc] - model_height_surface[iloc]);
      // Correct to observation height
      hofx[iloc] = model_T_surface[iloc] +
                   ufo::Constants::Lclr * (model_height_surface[iloc] - obs_height[iloc]);
    }
  }

  oops::Log::trace() << "airTemperatureAt2M_WRFDA::simobs done" << std::endl;
}

// ----------------------------------------
// Temperature operator using UKMO method
// ----------------------------------------

airTemperatureAt2M_UKMO::airTemperatureAt2M_UKMO(const std::string & name,
                                                 const Parameters_ & params)
  : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_temperature"));
  vars.push_back(oops::Variable("air_pressure"));
  vars.push_back(oops::Variable("air_pressure_at_surface"));
  requiredVars_ += vars;
}

void airTemperatureAt2M_UKMO::simobs(const ufo::GeoVaLs & gv,
                                     const ioda::ObsSpace & obsdb,
                                     std::vector<float> & hofx) const {
  oops::Log::trace() << "airTemperatureAt2M_UKMO::simobs start" << std::endl;

  // Create oops::Variable needed
  const oops::Variable model_height_var = oops::Variable(params_.geovarGeomZ.value());
  const oops::Variable model_p_var = oops::Variable("air_pressure");
  const oops::Variable model_p_surface_var = oops::Variable("air_pressure_at_surface");
  const oops::Variable model_T_var = oops::Variable("air_temperature");

  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const float missingFloat = util::missingValue<float>();
  const double height_used = 2000.0;
  const int model_nlevs = gv.nlevs(model_p_var);
  const double T_exponent = ufo::Constants::rd * ufo::Constants::Lclr / ufo::Constants::grav;

  // Create arrays needed
  std::vector<float> model_height_surface(nobs), model_p_surface(nobs),
      obs_height(nobs), lats(nobs);
  std::vector<double> model_p_2000m(nobs), model_T_2000m(nobs);
  obsdb.get_db("MetaData", "latitude", lats);
  bool convertLevel1GeopotentialHeight = false;

  // Get level 1 height.  If geopotential then convert to geometric height.
  if (isGeopotentialHeightVariable(params_.geovarGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs do geopotential conversion profile"
                       << std::endl;
    convertLevel1GeopotentialHeight = true;
  }

  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(model_height_surface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    oops::Log::trace()  << "ObsSfcCorrected::simulateObs do geopotential conversion surface"
                       << std::endl;
    convertToGeometricHeight(lats, model_height_surface);
  }

  // Read data in
  gv.get(model_p_surface, model_p_surface_var);
  obsdb.get_db("MetaData", params_.obsHeightName.value(), obs_height);

  // Loop to calculate hofx
  double model_T_surface;
  std::vector<double> profile_height(model_nlevs);
  std::vector<double> profile_pressure(model_nlevs);
  std::vector<double> profile_T(model_nlevs);
  int index = 0;
  double weight = 0.0;
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = missingFloat;
    if (obs_height[iloc] != missingFloat && model_p_surface[iloc] != missingFloat &&
        model_height_surface[iloc] != missingFloat) {
      // Get model data at this location
      gv.getAtLocation(profile_height, model_height_var, iloc);
      gv.getAtLocation(profile_pressure, model_p_var, iloc);
      gv.getAtLocation(profile_T, model_T_var, iloc);
      // Convert geopotential to geometric height if needed
      if (convertLevel1GeopotentialHeight) {
        convertToGeometricHeight(lats[iloc], profile_height);
      }
      vert_interp_weights_f90(model_nlevs, height_used, profile_height.data(), index, weight);

      // Vertical interpolation to get model pressure and temperature at 2000 m
      vert_interp_apply_f90(model_nlevs, profile_pressure.data(),
                            model_p_2000m[iloc], index, weight);
      vert_interp_apply_f90(model_nlevs, profile_T.data(),
                            model_T_2000m[iloc], index, weight);
      // Find model surface T using lowest T in model temperature profile
      model_T_surface = model_T_2000m[iloc] *
        std::pow((model_p_surface[iloc] / model_p_2000m[iloc]), T_exponent);
      // Correct to observation height
      hofx[iloc] = model_T_surface +
        ufo::Constants::Lclr * (model_height_surface[iloc] - obs_height[iloc]);
    }
  }
  oops::Log::trace() << "airTemperatureAt2M_UKMO::simobs done" << std::endl;
}

// ----------------------------------------
// Temperature operator using GSL method
// ----------------------------------------

airTemperatureAt2M_GSL::airTemperatureAt2M_GSL(const std::string & name,
                                               const Parameters_ & params)
  : SurfaceOperatorBase(name, params)
{
  oops::Variables vars;
  vars.push_back(oops::Variable(params_.geovarGeomZ.value()));
  vars.push_back(oops::Variable(params_.geovarSfcGeomZ.value()));
  vars.push_back(oops::Variable("air_temperature"));
  vars.push_back(oops::Variable("air_temperature_at_2m"));
  requiredVars_ += vars;
}

void airTemperatureAt2M_GSL::simobs(const ufo::GeoVaLs & gv,
                                    const ioda::ObsSpace & obsdb,
                                    std::vector<float> & hofx) const {
  oops::Log::trace() << "airTemperatureAt2M_GSL::simobs starting" << std::endl;
  // Setup parameters used throughout
  const size_t nobs = obsdb.nlocs();
  const float missingFloat = util::missingValue<float>();

  // Create arrays needed
  std::vector<float> lapse_rate(nobs), obs_height(nobs), lats(nobs),
                     model_height_surface(nobs), model_T_surface(nobs);
  obsdb.get_db("MetaData", "latitude", lats);

  // Check if GSL parameters exist
  if (!params_.gslParams.value()) {
    throw eckit::UserError("GSL correction requires 'gsl parameters' key to be provided", Here());
  }
  const auto& gsl_params = params_.gslParams.value().get();

  switch (gsl_params.temperatureLapseRateOption.value()) {
    case GslLapseRateOption::Constant: {
      const float lapse_rate_value =
          static_cast<float>(gsl_params.temperatureLapseRateValue.value() /
                              ufo::Constants::thousand);  // Convert K/km to K/m
      lapse_rate = std::vector<float>(nobs, lapse_rate_value);
      break;
    }
    case GslLapseRateOption::NoAdjustment: {
      lapse_rate = std::vector<float>(nobs, 0.0f);
      break;
    }
    case GslLapseRateOption::Local: {
      const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
      const int surface_level_index = gv.nlevs(geomz_var) - 1;
      lapse_rate = computeGslLocalLapseRate(gv, lats, params_.geovarGeomZ.value(), gsl_params,
                                            surface_level_index,
                                            "airTemperatureAt2M_GSL::simobs").lapseRate;
      break;
    }
  }

  // Get surface height.  If geopotential then convert to geometric height.
  gv.get(model_height_surface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
      oops::Log::trace() << "ObsSfcCorrected::simulateObs do geopotential conversion surface"
             << std::endl;
      convertToGeometricHeight(lats, model_height_surface);
  }

  // Read other data in
  if (gv.nlevs(oops::Variable{"air_temperature_at_2m"}) <= 0) {
     throw eckit::UserError("airTemperatureAt2M_GSL::simobs Variable "
                            "air_temperature_at_2m not present in Geoval", Here());
  }
  gv.get(model_T_surface, oops::Variable("air_temperature_at_2m"));
  obsdb.get_db("MetaData", params_.obsHeightName.value(), obs_height);

  // Loop to calculate hofx
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    hofx[iloc] = missingFloat;
    if (obs_height[iloc] != missingFloat && model_height_surface[iloc] != missingFloat &&
        lapse_rate[iloc] != missingFloat && model_T_surface[iloc] != missingFloat) {
      // Correct to observation height
      hofx[iloc] = model_T_surface[iloc] +
                   lapse_rate[iloc] * (model_height_surface[iloc] - obs_height[iloc]);
    }
  }

  oops::Log::trace() << "airTemperatureAt2M_GSL::simobs complete" << std::endl;
}

// -----------------------------------------------------------------------------

void airTemperatureAt2M_GSL::settraj(const GeoVaLs & gv,
                                     const ioda::ObsSpace & obsdb) {
  oops::Log::trace() << "airTemperatureAt2M_GSL::settraj start" << std::endl;

  if (!params_.gslParams.value()) {
    throw eckit::UserError("airTemperatureAt2M_GSL_TLAD requires 'gsl parameters'", Here());
  }
  const auto & gsl = params_.gslParams.value().get();
  needsLevelTemperature_ =
      (gsl.temperatureLapseRateOption.value() == GslLapseRateOption::Local);
  // Exact derivative of hofx = T_2m + lapse_rate(state) * (z_sfc - z_obs),
  // lapse_rate has zero derivative w.r.t. the model state in every case
  // except "Local" with "MatchSimobs" linearization
  // For Constant, NoAdjustment the dhofx reduces to dT_2m alone
  tlad_local_path_ = needsLevelTemperature_ &&
      (gsl.tladLapseRateOption.value() != GslTladLapseRateOption::NoAdjustment);
  // Optional: in the zero-derivative cases, also add dT_lev1 so the
  // increment still reaches the 3D field if air_temperature_at_2m carries
  // no independent active increment (see "tlad spread to lowest level")
  tlad_spread_to_lowest_level_ = gsl.tladSpreadToLowestLevel.value() && !tlad_local_path_;

  const size_t nobs = obsdb.nlocs();
  const float fmissing = util::missingValue<float>();

  std::vector<float> obs_height(nobs);
  std::vector<float> lats(nobs);
  std::vector<float> model_height_surface(nobs);
  obsdb.get_db("MetaData", "latitude", lats);
  obsdb.get_db("MetaData", params_.obsHeightName.value(), obs_height);

  // Surface model height (geopotential -> geometric if needed)
  gv.get(model_height_surface, oops::Variable(params_.geovarSfcGeomZ.value()));
  if (isGeopotentialHeightVariable(params_.geovarSfcGeomZ.value())) {
    convertToGeometricHeight(lats, model_height_surface);
  }

  // Per-location lapse rate, replicating simobs choices
  std::vector<double> lapse_rate(nobs, ufo::Constants::zero);
  std::vector<bool> lapse_rate_valid(nobs, true);
  lapse_rate_clamped_.assign(nobs, false);  // true where threshold clamp was active
  inv_lapse_layer_depth_.clear();
  lapse_top_level_index_ = -1;

  // surface_level_index_ is always needed: the Local lapse rate uses it, and
  // "tlad spread to lowest level" adds air_temperature at this level in TL/AD.
  const oops::Variable geomz_var = oops::Variable(params_.geovarGeomZ.value());
  surface_level_index_ = gv.nlevs(geomz_var) - 1;

  switch (gsl.temperatureLapseRateOption.value()) {
    case GslLapseRateOption::Constant: {
      const double lapse_rate_value =
          gsl.temperatureLapseRateValue.value() / ufo::Constants::thousand;  // K/km -> K/m
      std::fill(lapse_rate.begin(), lapse_rate.end(), lapse_rate_value);
      break;
    }
    case GslLapseRateOption::NoAdjustment: {
      std::fill(lapse_rate.begin(), lapse_rate.end(), ufo::Constants::zero);
      break;
    }
    case GslLapseRateOption::Local: {
      const GslLocalLapseRate local =
          computeGslLocalLapseRate(gv, lats, params_.geovarGeomZ.value(), gsl,
                                   surface_level_index_, "airTemperatureAt2M_GSL::settraj");
      lapse_top_level_index_ = local.lapseTopLevelIndex;
      lapse_rate_clamped_ = local.lapseRateClamped;

      if (needsLevelTemperature_) {
        inv_lapse_layer_depth_.assign(nobs, ufo::Constants::zero);
      }
      for (size_t iobs = 0; iobs < nobs; ++iobs) {
        if (local.lapseRate[iobs] == fmissing) {
          lapse_rate_valid[iobs] = false;
          continue;
        }
        lapse_rate[iobs] = static_cast<double>(local.lapseRate[iobs]);
        if (needsLevelTemperature_) {
          inv_lapse_layer_depth_[iobs] =
              ufo::Constants::one / static_cast<double>(local.lapseLayerDepth[iobs]);
        }
      }
      break;
    }
  }

  // Save trajectory state.
  dz_sfc_minus_obs_.assign(nobs, ufo::Constants::zero);
  valid_.assign(nobs, false);
  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    if (obs_height[iloc] != fmissing && model_height_surface[iloc] != fmissing &&
        lapse_rate_valid[iloc]) {
      dz_sfc_minus_obs_[iloc] =
          static_cast<double>(model_height_surface[iloc]) - static_cast<double>(obs_height[iloc]);
      valid_[iloc] = true;
    }
  }

  oops::Log::trace() << "airTemperatureAt2M_GSL::settraj done" << std::endl;
}

// -----------------------------------------------------------------------------

void airTemperatureAt2M_GSL::getLevelIncrements(const GeoVaLs & gv,
                                                std::vector<double> & dT_lev1,
                                                std::vector<double> & dT_top) const {
  const size_t nobs = valid_.size();
  const oops::Variable model_T_var = oops::Variable("air_temperature");
  if (tlad_local_path_) {
    dT_lev1.assign(nobs, ufo::Constants::zero);
    dT_top.assign(nobs, ufo::Constants::zero);
    gv.getAtLevel(dT_lev1, model_T_var, surface_level_index_);
    gv.getAtLevel(dT_top, model_T_var, lapse_top_level_index_);
  } else if (tlad_spread_to_lowest_level_) {
    dT_lev1.assign(nobs, ufo::Constants::zero);
    gv.getAtLevel(dT_lev1, model_T_var, surface_level_index_);
  }
}

// -----------------------------------------------------------------------------

void airTemperatureAt2M_GSL::TL(const GeoVaLs & gv,
                                std::vector<double> & dhofx) const {
  oops::Log::trace() << "airTemperatureAt2M_GSL::TL start" << std::endl;
  const size_t nobs = dhofx.size();
  const double dmissing = util::missingValue<double>();
  const oops::Variable T2m_var = oops::Variable("air_temperature_at_2m");

  std::vector<double> dT_2m(nobs);
  gv.get(dT_2m, T2m_var);

  std::vector<double> dT_lev1;
  std::vector<double> dT_top;
  getLevelIncrements(gv, dT_lev1, dT_top);

  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    if (!valid_[iloc]) { dhofx[iloc] = dmissing; continue; }

    double v = dT_2m[iloc];
    if (tlad_local_path_ && !lapse_rate_clamped_[iloc]) {
      const double dlapse = (dT_lev1[iloc] - dT_top[iloc]) * inv_lapse_layer_depth_[iloc];
      v += dlapse * dz_sfc_minus_obs_[iloc];
    } else if (tlad_spread_to_lowest_level_) {
      v += dT_lev1[iloc];
    }
    dhofx[iloc] = v;
  }

  oops::Log::trace() << "airTemperatureAt2M_GSL::TL done" << std::endl;
}

// -----------------------------------------------------------------------------

void airTemperatureAt2M_GSL::AD(GeoVaLs & gv,
                                std::vector<double> & dhofx) const {
  oops::Log::trace() << "airTemperatureAt2M_GSL::AD start" << std::endl;
  const size_t nobs = dhofx.size();
  const double dmissing = util::missingValue<double>();
  const oops::Variable model_T_var = oops::Variable("air_temperature");
  const oops::Variable T2m_var = oops::Variable("air_temperature_at_2m");

  std::vector<double> dT_2m(nobs);
  gv.get(dT_2m, T2m_var);

  std::vector<double> dT_lev1;
  std::vector<double> dT_top;
  getLevelIncrements(gv, dT_lev1, dT_top);

  for (size_t iloc = 0; iloc < nobs; ++iloc) {
    if (!valid_[iloc] || dhofx[iloc] == dmissing) continue;
    const double dy = dhofx[iloc];
    dT_2m[iloc] += dy;
    if (tlad_local_path_ && !lapse_rate_clamped_[iloc]) {
      const double w = dy * dz_sfc_minus_obs_[iloc] * inv_lapse_layer_depth_[iloc];
      dT_lev1[iloc] += w;
      dT_top[iloc]  -= w;
    } else if (tlad_spread_to_lowest_level_) {
      dT_lev1[iloc] += dy;
    }
  }

  gv.putAtLevel(dT_2m, T2m_var, 0);
  if (tlad_local_path_) {
    gv.putAtLevel(dT_lev1, model_T_var, surface_level_index_);
    gv.putAtLevel(dT_top, model_T_var, lapse_top_level_index_);
  } else if (tlad_spread_to_lowest_level_) {
    gv.putAtLevel(dT_lev1, model_T_var, surface_level_index_);
  }

  oops::Log::trace() << "airTemperatureAt2M_GSL::AD done" << std::endl;
}

}  // namespace ufo
