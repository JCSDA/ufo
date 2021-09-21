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
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorSfcPressure> makerSteps_("ObsErrorFactorSfcPressure");

// -----------------------------------------------------------------------------

ObsErrorFactorSfcPressure::ObsErrorFactorSfcPressure(const eckit::Configuration &config)
  : invars_() {
  oops::Log::debug() << "ObsErrorFactorSfcPressure: config = " << config << std::endl;
  const float tiny_float = FLT_MIN;
  const float huge_float = FLT_MAX;
  const float missing = util::missingValue(missing);
  // Initialize options
  options_.reset(new ObsErrorFactorSfcPressureParameters());
  options_->deserialize(config);

  // Initialize error_min, max from options. Make sure they are sane.
  const float error_min = options_->error_min.value();
  const float error_max = options_->error_max.value();
  ASSERT(error_min < error_max);

  // The starting (un-inflated) value of obserror. If running in sequence of filters,
  // then it is probably found in ObsErrorData, otherwise, it is probably ObsError.
  const std::string errgrp = options_->original_obserr.value();
  invars_ += Variable("surface_pressure@"+errgrp);

  // Include list of required data from ObsValue
  invars_ += Variable("surface_pressure@ObsValue");
  invars_ += Variable("air_temperature@ObsValue");

  // Include list of required data from MetaData
  invars_ += Variable("station_elevation@MetaData");

  // Include list of required data from GeoVaLs
  invars_ += Variable("surface_pressure@GeoVaLs");
  invars_ += Variable("air_pressure@GeoVaLs");
  const std::string geovar_temp = options_->geovar_temp.value();
  invars_ += Variable(geovar_temp + "@GeoVaLs");
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  invars_ += Variable(geovar_sfc_geomz + "@GeoVaLs");
}

// -----------------------------------------------------------------------------

ObsErrorFactorSfcPressure::~ObsErrorFactorSfcPressure() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorSfcPressure::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue(missing);
  static constexpr float g_over_rd = 1.0f*Constants::grav/Constants::rd;
  const float lapse_rate = 1.0f*Constants::Lclr;

  // Ensure that only one output variable is expected.
  ASSERT(obserr.nvars() == 1);

  // Get dimensions
  size_t nlocs = data.nlocs();
  size_t nlevs = data.nlevs(Variable("air_pressure@GeoVaLs"));

  // Get min, max error values
  const float error_min = options_->error_min.value();
  const float error_max = options_->error_max.value();

  // Get MetaData of station elevation
  std::vector<float> ob_elevation(nlocs);
  data.get(Variable("station_elevation@MetaData"), ob_elevation);

  // Get ObsValue of surface pressure and temperature (possibly missing).
  std::vector<float> ob_pressure_sfc(nlocs);
  data.get(Variable("surface_pressure@ObsValue"), ob_pressure_sfc);
  std::vector<float> ob_temperature_sfc(nlocs);
  data.get(Variable("air_temperature@ObsValue"), ob_temperature_sfc);

  // Get original ObsError of surface pressure
  std::vector<float> currentObserr(nlocs);
  const std::string errgrp = options_->original_obserr.value();
  data.get(Variable("surface_pressure@"+errgrp), currentObserr);

  // Get GeoVals of surface altitude, pressure, and temperature.
  std::vector<float> model_elevation(nlocs);
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  data.get(Variable(geovar_sfc_geomz + "@GeoVaLs"), model_elevation);
  std::vector<float> model_pres_sfc(nlocs);
  data.get(Variable("surface_pressure@GeoVaLs"), model_pres_sfc);

  // Get GeoVals of air pressure [Pa] and temperature in vertical column.
  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    data.get(Variable("air_pressure@GeoVaLs"), level, prsl[ilev]);
  }
  std::vector<std::vector<float>> tair(nlevs, std::vector<float>(nlocs));
  const std::string geovar_temp = options_->geovar_temp.value();
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    data.get(Variable(geovar_temp + "@GeoVaLs"), level, tair[ilev]);
  }

  // TODO(gthompsn): model temp is virtual_temp whereas obs is sensible temp.
  // Also, for now, assigning model lowest level temp as surface temp.
  // Lastly, need to make positive that assumption of 0th level is nearest sfc.
  std::vector<float> model_temp_sfc(nlocs);
  model_temp_sfc = tair[0];

  float obserror, new_error, error_factor;
  float rdelz, rdp, drdp, ddiff, tges, pges, tges2, pges2, psges, pgesorig, drbx;
  int iv = 0;

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // If missing station_elevation, set error factor to something very high (for rejection).
    if (ob_elevation[iloc] == missing) {
       obserr[iv][iloc] = 99.9f;
    } else {
      rdelz = ob_elevation[iloc]-model_elevation[iloc];

      // If more than a km between model and real elevation, set error factor linearly higher.
      if (std::abs(rdelz) > 5000.0f) {
       obserr[iv][iloc] = 50.0f;
      } else if (std::abs(rdelz) > 1000.0f) {
       obserr[iv][iloc] = 3.0f + 47.0f*(std::abs(rdelz)-1000.0f)/4000.0f;
      } else {
        pgesorig = model_pres_sfc[iloc]*0.001;             // Converting Pascals to cb
        psges = log(pgesorig);

        // Calculating drbx with observed temperature.
        // If ob temperature missing, then check if model ground is above or below actual ground.
        if (ob_temperature_sfc[iloc] != missing) {
          drbx = 0.5f*std::abs(model_temp_sfc[iloc]-ob_temperature_sfc[iloc])
                 +0.2f+0.005f*std::abs(rdelz);
          tges = 0.5f*(model_temp_sfc[iloc]+ob_temperature_sfc[iloc]);
        } else {
          // TODO(gthompsn): If model terrain below real world, grab nearest Temp,Pres from a
          // vertical profile. Shortcut for now is assume lapse_rate and hydrostatic assumption
          // over rdelz. The 2.5 addition is **arbitrary** and lapse rate assumption is 5C/km.
          if (std::abs(rdelz) < 5.0) {
            tges = model_temp_sfc[iloc];
            drbx = 0.1;
          } else if (rdelz > 0.0) {
            tges2 = model_temp_sfc[iloc] - lapse_rate*rdelz;
            drbx = 0.5f*std::abs(model_temp_sfc[iloc]-tges2)+2.5f+0.005f*std::abs(rdelz);
            tges = 0.5f*(model_temp_sfc[iloc]+tges2);
          } else {
            tges = model_temp_sfc[iloc] - 0.5f*lapse_rate*rdelz;
            tges2 = tges - lapse_rate*rdelz;
            drbx = 0.5f*std::abs(model_temp_sfc[iloc]-tges2)+2.5f+0.005f*std::abs(rdelz);
            drbx = drbx - 0.5f*lapse_rate*rdelz;
          }
        }

        rdp = g_over_rd*rdelz/tges;
        pges = exp(log(pgesorig) - rdp);
        drdp = pges*(g_over_rd*std::abs(rdelz)*drbx/(tges*tges));
        ddiff = ob_pressure_sfc[iloc]*0.001f - pges;        // innovation in cb

        // make adjustment to observational error (also convert to cb)
        obserror = currentObserr[iloc]*0.001f;
        // TODO(gthompsn): Maybe reduce obserror by 0.7 for data near sea-level and small delta-Z.
        // if (ob_elevation[iloc] < 10.0f && rdelz < 5.0f) {
        //   obserror = obserror*0.7;
        // }
        new_error = obserror + drdp;
        new_error = std::max(error_min*0.001f, std::min(new_error, error_max*0.001f));
        error_factor = std::max(0.7f, new_error/obserror);

        obserr[iv][iloc] = error_factor;
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSfcPressure::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
