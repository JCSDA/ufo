/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSfcPressure.h"

#include <float.h>

#include <algorithm>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorSfcPressure> makerSteps_("ObsErrorFactorSfcPressure");

// -----------------------------------------------------------------------------

ObsErrorFactorSfcPressure::ObsErrorFactorSfcPressure(const eckit::Configuration &config)
  : invars_() {
  oops::Log::debug() << "ObsErrorFactorSfcPressure: config = " << config << std::endl;
  const float tiny_float = FLT_MIN;
  const float huge_float = FLT_MAX;
  const float missing = util::missingValue<float>();
  // Initialize options
  options_.reset(new ObsErrorFactorSfcPressureParameters());
  options_->deserialize(config);

  // The starting (un-inflated) value of obserror. If running in sequence of filters,
  // then it is probably found in ObsErrorData, otherwise, it is probably ObsError.
  const std::string errgrp = options_->original_obserr.value();
  invars_ += Variable(errgrp+"/stationPressure");

  // Include list of required data from ObsValue
  invars_ += Variable("ObsValue/stationPressure");
  invars_ += Variable("ObsValue/virtualTemperature");
  invars_ += Variable("ObsValue/airTemperature");

  // Include list of required data from MetaData
  invars_ += Variable("MetaData/stationElevation");

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/surface_pressure");
  invars_ += Variable("GeoVaLs/air_pressure");
  invars_ += Variable("GeoVaLs/virtual_temperature");
  const std::string geovar_geomz = options_->geovar_geomz.value();
  invars_ += Variable("GeoVaLs/" + geovar_geomz);
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  invars_ += Variable("GeoVaLs/" + geovar_sfc_geomz);
}

// -----------------------------------------------------------------------------

ObsErrorFactorSfcPressure::~ObsErrorFactorSfcPressure() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorSfcPressure::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue<float>();
  static constexpr float g_over_rd = 1.0f*Constants::grav/Constants::rd;
  const float lapse_rate = 1.0f*Constants::Lclr;

  // If no observations on this processor then nothing to do
  if (data.nlocs() == 0) return;

  // Ensure that only one output variable is expected.
  ASSERT(obserr.nvars() == 1);

  // Get dimensions
  size_t nlocs = data.nlocs();
  size_t nlevs = data.nlevs(Variable("GeoVaLs/air_pressure"));

  // Get MetaData of station elevation
  std::vector<float> ob_elevation(nlocs);
  data.get(Variable("MetaData/height"), ob_elevation);

  // Get ObsValue of surface pressure
  std::vector<float> ob_pressure_sfc(nlocs);
  data.get(Variable("ObsValue/stationPressure"), ob_pressure_sfc);

  // Get ObsValue of virtual temperature (optional), initialize
  // the vector as missing value and get values only if it exists
  std::vector<float> ob_temp_sfc(nlocs, missing);
  if (data.has(Variable("ObsValue/virtualTemperature"))) {
    data.get(Variable("ObsValue/virtualTemperature"), ob_temp_sfc);
  }
  if (data.has(Variable("ObsValue/airTemperature"))) {
    std::vector<float> ob_temp_sfc2(nlocs);
    data.get(Variable("ObsValue/airTemperature"), ob_temp_sfc2);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (ob_temp_sfc[iloc] == missing &&
          ob_temp_sfc2[iloc] != missing) {
        ob_temp_sfc[iloc] =  ob_temp_sfc2[iloc];
      }
    }
  }

  // Get original ObsError of surface pressure
  std::vector<float> currentObserr(nlocs);
  const std::string errgrp = options_->original_obserr.value();
  data.get(Variable(errgrp + "/stationPressure"), currentObserr);

  // Get GeoVaLs of surface altitude and pressure
  std::vector<float> model_elevation(nlocs);
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  data.get(Variable("GeoVaLs/" + geovar_sfc_geomz), model_elevation);
  if (geovar_sfc_geomz.find("geopotential_height") != std::string::npos) {
    // Transform geopotential height to geometric height
    std::vector<float> latitude(nlocs);
    data.get(Variable("MetaData/latitude"), latitude);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      model_elevation[iloc] = formulas::Geopotential_to_Geometric_Height(latitude[iloc],
                              model_elevation[iloc]);
    }
  }
  std::vector<float> model_pres_sfc(nlocs);
  data.get(Variable("GeoVaLs/surface_pressure"), model_pres_sfc);

  // Get GeoVaLs pointer for retrieving vertical profiles
  const ufo::GeoVaLs * gvals = data.getGeoVaLs();
  std::vector<double> temperature_gval(nlevs, 0.0);
  std::vector<double> pressure_gval(nlevs, 0.0), logp(nlevs);

  // Get index of model bottom level
  int levbot = nlevs-1;
  gvals->getAtLocation(pressure_gval, oops::Variable{"air_pressure"}, 0);
  if (pressure_gval[0] > pressure_gval[nlevs-1]) levbot = 0;

  // Get GeoVaLs of altitude and virtual temperature at model bottom level
  std::vector<float> model_temp_sfc(nlocs);
  data.get(Variable("GeoVaLs/virtual_temperature"), levbot, model_temp_sfc);
  std::vector<float> model_elevation_bot(nlocs);
  const std::string geovar_geomz = options_->geovar_geomz.value();
  data.get(Variable("GeoVaLs/" + geovar_geomz), levbot, model_elevation_bot);
  if (geovar_geomz.find("geopotential_height") != std::string::npos) {
    // Transform geopotential height to geometric height
    std::vector<float> latitude(nlocs);
    data.get(Variable("MetaData/latitude"), latitude);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      model_elevation_bot[iloc] = formulas::Geopotential_to_Geometric_Height(
                   latitude[iloc], model_elevation_bot[iloc]);
    }
  }

  // Derive error inflation based on surface pressure correction and
  // further inflate error by sqrt(factor_dup)
  float new_error;
  float rdelz, drdp, tges, tges2, pges, drbx, logp_sfc;
  int iv = 0;

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // If missing surface pressure or station_elevation,
    // set error factor to something very high (for rejection).
    if (ob_elevation[iloc] == missing || ob_pressure_sfc[iloc] == missing) {
      obserr[iv][iloc] = 1.0f;
    } else {
      rdelz = ob_elevation[iloc]-model_elevation[iloc];
      // Calculate drbx with observed surface temperature (if not missing)
      // OR model temperature at observed surface altitude
      if ((ob_temp_sfc[iloc] < 350.0f) && (ob_temp_sfc[iloc] > 150.0f)) {
        // drbx is the amount of temperature changes from temp diff & altitude diff
        drbx = 0.5f*std::abs(model_temp_sfc[iloc]-ob_temp_sfc[iloc])
               +0.2f+0.005f*std::abs(rdelz);
        tges = 0.5f*(model_temp_sfc[iloc]+ob_temp_sfc[iloc]);
      } else {
        // Get model temperature at obs_elevation (tges2)
        // Interpolate GeoVaLs of virtual_temperature to obs_elevation
        gvals->getAtLocation(pressure_gval, oops::Variable{"air_pressure"}, iloc);
        gvals->getAtLocation(temperature_gval, oops::Variable{"virtual_temperature"}, iloc);
        for (size_t ilev = 0; ilev < nlevs; ++ilev) logp[ilev] = std::log(pressure_gval[ilev]);
        ufo::PiecewiseLinearInterpolation vert_interp_model(logp, temperature_gval);
        logp_sfc = std::log(ob_pressure_sfc[iloc]);
        // Adjust drbx & tges
        tges2 = vert_interp_model(logp_sfc);
        tges = 0.5f*(model_temp_sfc[iloc]+tges2);
        drbx = 0.5f*std::abs(model_temp_sfc[iloc]-tges2)
               +2.5f+0.005f*std::abs(rdelz);
        if (rdelz < 0.0f) {
          // Extrapolate GeoVaLs of surface virtual_temperature to obs_elevation
          tges = tges - (lapse_rate * 0.5f * rdelz);
          drbx = drbx - (lapse_rate * 0.5f * rdelz);
        }
      }
      // Adjust model surface pressure
      pges = std::exp(std::log(model_pres_sfc[iloc]) - (rdelz * g_over_rd / tges));
      // Adjust observation error
      drdp = pges * (g_over_rd * std::abs(rdelz) * drbx / (tges * tges));
      new_error = (currentObserr[iloc] + drdp);
      obserr[iv][iloc] = new_error / currentObserr[iloc];
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSfcPressure::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
