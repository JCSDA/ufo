/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/StableLayersCloudTopPressure.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/ufo_utils.interface.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<StableLayersCloudTopPressure>
 makerStableLayersCloudTopPressure_("StableLayersCloudTopPressure");

// -----------------------------------------------------------------------------

StableLayersCloudTopPressure::StableLayersCloudTopPressure
(const eckit::LocalConfiguration & conf)
  : invars_(), channels_() {
  oops::Log::trace() << "StableLayersCloudTopPressure constructor" << std::endl;
  // Validate and deserialize options.
  options_.validateAndDeserialize(conf);

  // Get channels used
  const std::set<int> chanset = oops::parseIntSet(options_.channels.value());
  channels_.assign(chanset.begin(), chanset.end());
  ASSERT(channels_.size() > 0);

  // Obtain variables associated with the `where` options.
  invars_ += getAllWhereVariables(options_.where);

  // Required tropopause level.
  invars_ += Variable(options_.tropopauseLevel.value());

  // Required bias corrected brightness temperature.
  invars_ += Variable(options_.correctedBrightnessTemperature.value());

  // GeoVaLs.
  invars_ += Variable(std::string("GeoVaLs/air_temperature"));
  invars_ += Variable(std::string("GeoVaLs/air_pressure"));
  invars_ += Variable(std::string("GeoVaLs/relative_humidity"));

  // Model brightness temperature.
  invars_ += Variable("ObsDiag/brightness_temperature_from_atmosphere_layer_to_toa", channels_);
}

// -----------------------------------------------------------------------------

StableLayersCloudTopPressure::~StableLayersCloudTopPressure() {
  oops::Log::trace() << "StableLayersCloudTopPressure destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void StableLayersCloudTopPressure::compute(const ObsFilterData & in,
                                        ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "StableLayersCloudTopPressure compute start" << std::endl;
  // Missing float value.
  const float missing = util::missingValue<float>();

  // ObsSpace.
  ioda::ObsSpace & obsdb = in.obsspace();

  // Number of locations.
  const size_t nlocs = in.nlocs();

  // Number of channels.
  const size_t nchans = obsdb.nchans();

  // GeoVaLs.
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Number of model levels.
  const int nlevs = in.nlevs(Variable("GeoVaLs/air_pressure"));

  // Flag to use brightness temperature constraint method.
  const bool useBTConstraint = options_.useBTConstraint.value();

  // Constants.
  const float tempLimitWarm = options_.tempLimitWarm.value();
  const float tempLimitCold = options_.tempLimitCold.value();
  const float stableDensity = options_.stableDensity.value();
  const float relativeHumidityDensity = options_.relativeHumidityDensity.value();
  const double relativeHumidityOffset = options_.relativeHumidityOffset.value();
  const double relativeHumidityMinimum = options_.relativeHumidityMinimum.value();

  // Channel index.
  const int windowChannel = options_.windowChannel.value();

  // Tropopause level.
  std::vector<int> tropopauseLevel(nlocs);

  // Bias corrected brightness temperature.
  std::vector<float> correctedBrightnessTemperature(nlocs);

  // Get required input variables.
  in.get(Variable(options_.tropopauseLevel.value()), tropopauseLevel);
  in.get(Variable(options_.correctedBrightnessTemperature.value()), correctedBrightnessTemperature);

  // Vectors of GeoVaLs.
  std::vector<double> gv_p(nlevs);
  std::vector<double> gv_t(nlevs);
  std::vector<double> gv_rh(nlevs);

  // Get the overcast brightness temperature for the window channel at every level and for
  // every location.
  const std::string btOvercastName =
                  std::string("ObsDiag/brightness_temperature_from_atmosphere_layer_to_toa")
                                      + "_" + std::to_string(windowChannel);
  std::vector<std::vector<float>> overcast_bt(nlevs, std::vector<float>(nlocs));
  for (int ilev = 0; ilev < nlevs; ++ilev) {
    in.get(Variable(btOvercastName), ilev, overcast_bt[ilev]);
  }

  // Assign vectors to write flags/output to.
  std::vector<bool> nostablelayers_flag;
  std::vector<float> bt_cloud_top(nlocs);
  std::vector<float> cloud_top_pressure(nlocs);
  std::vector<float> standard_deviation(nlocs);

  // Get the overcast brightness temperature and standard deviation from the ObsSpace.
  // If this is the first call to the obsfunction, the vector will be filled with missing
  // values. If it is not the first call, the vector will be filled with the values from the
  // previous call. This is to ensure that observations which have already been assigned a
  // brightness temperature and standard deviation that fulfil the requirements to be a stable
  // layer are not overwritten. Note that if the overcast brightness temperature is present,
  // the standard deviation will also be present, so there is no need to test for both.
  if (obsdb.has("MetaData", "overcastBriTempAtCloudTop")) {
    obsdb.get_db("MetaData", "overcastBriTempAtCloudTop", bt_cloud_top);
    obsdb.get_db("MetaData", "standardDeviationAtCloudTop", standard_deviation);
  } else {
    bt_cloud_top.assign(nlocs, missing);
    standard_deviation.assign(nlocs, missing);
  }

  // Get the nostablelayers flag from the ObsSpace. If this is the first call to the obsfunction,
  // the values must have been created with a Create Diagnostic Flags filter. If it is not the
  // first call, the values will be filled with the values from the previous call.
  // nostablelayers_flag has dimension (nlocs*nchans).
  obsdb.get_db("DiagnosticFlags/NoStableLayers", "brightnessTemperature", nostablelayers_flag);

  // Vector of locations that pass the 'where' clause in the sample
  // (all true if there is no where clause).
  const std::vector<bool> apply = processWhere(options_.where, in, options_.whereOperator);

  // Loop over locations.
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Skip this location if it satisfies the requirements to be a stable layer.
    // This is set according to the contents of the where block. In a typical use case, this will
    // happen if it has already been assigned a brightness temperature in a previous call and if
    // the brightness temperature constraint is set to false.
    // The location will keep its previous value of cloud top pressure, assigned via out[0],
    // as well as brightness temperature and nostablelayers_flag, assigned above.
    // This is to avoid overwriting the brightness temperature at cloud top with missing values.
    if (!apply[jloc]) {
      continue;
    }

    // We also skip this location if the bias corrected brightness temperature is missing. In this
    // case we still need to output the cloud_top_pressure via out[0] so that the ObsSpace receives
    // a value of missing. The brightness temperature at cloud top will be output as missing after
    // the loop. Since a missing value of bias corrected brightness temperature prevents us from
    // searching for stable layers, we set nostablelayers_flag to true for this location.
    if (correctedBrightnessTemperature[jloc] == missing) {
      // Set nostablelayers_flag to true for the entire row (every channel) corresponding to jloc.
      for (size_t i = jloc * nchans; i < (jloc + 1) * nchans; ++i) {
        nostablelayers_flag[i] = true;
      }
      // Set the cloud top pressure at this location to missing.
      out[0][jloc] = missing;
      continue;
    }

    // If the tropopause level is the top of the atmosphere, reset it so that it is one level
    // below the top. This is to avoid out-of-bounds access in the subsequent calculations.
    // The cloud top pressure is very unlikely to be at the top of the atmosphere, so this is a
    // reasonable workaround.
    if (tropopauseLevel[jloc] == 0) {
      tropopauseLevel[jloc] = 1;
    }

    // Vector to populate with the stability weights.
    std::vector<double> stability_weights(nlevs, 0.0);

    // Get GeoVaLs at this location.
    gv->getAtLocation(gv_p, oops::Variable{"air_pressure"}, jloc);
    gv->getAtLocation(gv_t, oops::Variable{"air_temperature"}, jloc);
    gv->getAtLocation(gv_rh, oops::Variable{"relative_humidity"}, jloc);

    // Log of model pressure.
    std::vector<double> gv_log_p(nlevs);
    std::transform(gv_p.cbegin(), gv_p.cend(), gv_log_p.begin(),
                  [](double p) -> double {return std::log(p);});

    // Proceed through the atmosphere from the bottom up to the tropopause.
    for (int ilev = nlevs-1; ilev >= tropopauseLevel[jloc]; --ilev) {
      // -----------------------------------------------------------------------------
      // In this section, we test for the presence of stable layers.

      // Calculate the lapse rate dT/dlogp for this layer.
      const double lapse_rate = (gv_t[ilev-1] - gv_t[ilev]) /
                                          (gv_log_p[ilev-1] - gv_log_p[ilev]);

      // Calculate the mean saturated adiabatic lapse rate for this layer and convert to
      // log(p) coordinates.
      const double sat_lapse_rate_middle = computeSALR(gv_p[ilev], gv_t[ilev]);
      const double sat_lapse_rate_upper = computeSALR(gv_p[ilev-1], gv_t[ilev-1]);
      const double sat_lapse_rate = (gv_p[ilev]*sat_lapse_rate_middle +
                                              gv_p[ilev-1]*sat_lapse_rate_upper) / 2.0;

      // Conduct the three stable layers tests.
      // Test 1: Check for layer static stability. If not absolutely stable or saturated neutral,
      //         then first test is failed.
      if (lapse_rate > sat_lapse_rate) {
        continue;
      }

      // The next two tests are only conducted if the brightness temperature constraint flag is
      // set on input.
      if (useBTConstraint) {
        const double bt_diff = correctedBrightnessTemperature[jloc] - overcast_bt[ilev][jloc];
        // Test 2: Compare bias corrected brightness temperature with level brightness temperature.
        //         If the former is too much warmer than the latter, the test is failed.
        if (bt_diff > tempLimitWarm) {
          continue;
        }

        // Test 3: Compare bias corrected brightness temperature with level brightness temperature.
        //         If the former is too much colder than the latter, the test is failed.
        if (bt_diff < tempLimitCold) {
          continue;
        }
      }

      // -----------------------------------------------------------------------------
      // In this section, we calculate the stability weights for each level if all tests have
      // passed. The way this is done depends on both the level and the brightness temperature
      // constraint flag.
      // If we are not at the bottom of the atmosphere, the weight assigned to this level depends
      // on how stable the layer below is.
      // If the brightness temperature constraint flag is set, the stability weighting also takes
      // account of the brightness temperature difference.

      // Calculate layer temperature mean for stability weighting.
      const double mean_temp = (gv_t[ilev] + gv_t[ilev-1]) / 2.0;

      // Calculate the weightings that correspond to the lapse rate and relative humidity.
      const double lapse_weight = std::max(std::min((sat_lapse_rate - lapse_rate) /
                                                    (stableDensity * mean_temp), 1.0), 0.0);
      const double rh_weight = std::max(std::min((gv_rh[ilev] / relativeHumidityDensity) +
                                  relativeHumidityOffset, 1.0), relativeHumidityMinimum);

      // Initialise the lapse weight for the layer below to 0.0 and calculate its value if we are
      // not at the bottom of the atmosphere.
      double lapse_weight_below = 0.0;
      if (ilev < nlevs-1) {
        // Calculate layer temperature and pressure mean for the layer below.
        double mean_temp_below = (gv_t[ilev+1] + gv_t[ilev]) / 2.0;

        // Calculate the lapse rate dT/dlogp for the layer below.
        double lapse_rate_below = (gv_t[ilev] - gv_t[ilev+1]) / (gv_log_p[ilev] - gv_log_p[ilev+1]);

        // Calculate the mean saturated adiabatic lapse rate for the layer below and convert
        // to log(p) coordinates.
        const double sat_lapse_rate_lower = computeSALR(gv_p[ilev+1], gv_t[ilev+1]);
        const double sat_lapse_rate_below = (gv_p[ilev+1]*sat_lapse_rate_lower +
                                            gv_p[ilev]*sat_lapse_rate_middle) / 2.0;

        // Calculate the lapse weight for the layer below.
        lapse_weight_below = std::max(std::min((sat_lapse_rate_below - lapse_rate_below) /
                                        (stableDensity * mean_temp_below), 1.0), 0.0);
      }

      if (!useBTConstraint) {
        // Calculate the stability weight for each level. If we are not at the bottom of
        // the atmosphere, the weight assigned to this level depends on how stable the layer
        // below is. The current level will be penalised if the layer below is also very stable.
        if (ilev < nlevs-1) {
          stability_weights[ilev] = lapse_weight * (1.0 - lapse_weight_below) * rh_weight;
        } else {
          stability_weights[ilev] = lapse_weight * rh_weight;
        }
      } else {
        // If the brightness temperature constraint flag is set, the stability weighting also
        // takes account of the brightness temperature difference.
        const double bt_diff = correctedBrightnessTemperature[jloc] - overcast_bt[ilev][jloc];
        const double temp_weight = bt_diff < 0.0 ?
                                1.0 - bt_diff / tempLimitCold : 1.0 - bt_diff / tempLimitWarm;

        if (ilev < nlevs-1) {
          stability_weights[ilev] = lapse_weight * (1.0 - lapse_weight_below) * rh_weight
                                                                                * temp_weight;
        } else {
          stability_weights[ilev] = lapse_weight * rh_weight * temp_weight;
        }
      }
    }

    // -----------------------------------------------------------------------------
    // In this section, we fit a quadratic curve to calculate the true cloud top pressure.

    // Find the index of the maximum element in stability_weights. We use the reverse iterator
    // to find the last occurrence of the maximum element. This is because we want to find the
    // maximum that corresponds to the lowest level in the atmosphere. The stable layers scheme
    // is generally used for low-cloud situations, so given a choice between two levels with the
    // same stability weight, we choose the lower level.
    const auto max_weight = std::max_element(stability_weights.rbegin(), stability_weights.rend());
    const int max_weight_index = std::distance(stability_weights.begin(), max_weight.base()) - 1;

    // Set nostablelayers_flag to true at this location if the maximum stability weight is
    // below the threshold.
    if ((!useBTConstraint && *max_weight < 0.2) ||
        (useBTConstraint && *max_weight < 0.01)) {
      // Set nostablelayers_flag to true for the entire row (every channel) corresponding to jloc.
      for (size_t i = jloc * nchans; i < (jloc + 1) * nchans; ++i) {
        nostablelayers_flag[i] = true;
      }
      // Set both the cloud top pressure and the brightness temperature at this location to missing.
      out[0][jloc] = missing;
      bt_cloud_top[jloc] = missing;
      standard_deviation[jloc] = missing;
      continue;
    }

    // Calculate the cloud top pressure. If we are at the bottom of the atmosphere or at the
    // level of the tropopause, the cloud top pressure is assigned as the pressure at this level.
    // If we are at any other level, we fit a quadratic to find the true maximum.
    if (max_weight_index == (nlevs -1) || max_weight_index == tropopauseLevel[jloc]) {
      cloud_top_pressure[jloc] = gv_p[max_weight_index];
    } else {
      // Find the true cloud top pressure by fitting a quadratic to the level of the maximum
      // stability weight and its two neighbours, and then finding the maximum of the curve.
      const double log_pressure_diff_upper = gv_log_p[max_weight_index - 1] -
                                                                    gv_log_p[max_weight_index];
      const double log_pressure_diff_lower = gv_log_p[max_weight_index + 1] -
                                                                    gv_log_p[max_weight_index];
      const double weight_diff_upper = stability_weights[max_weight_index - 1] -
                                                           stability_weights[max_weight_index];
      const double weight_diff_lower = stability_weights[max_weight_index + 1] -
                                                           stability_weights[max_weight_index];
      const double numerator = std::pow(log_pressure_diff_lower, 2) * weight_diff_upper -
                                      std::pow(log_pressure_diff_upper, 2) * weight_diff_lower;
      const double denominator = 2 * (log_pressure_diff_lower * weight_diff_upper -
                                                  log_pressure_diff_upper * weight_diff_lower);
      if (denominator == 0.0) {
        cloud_top_pressure[jloc] = gv_p[max_weight_index];
      } else {
        cloud_top_pressure[jloc] = gv_p[max_weight_index] * std::exp(numerator / denominator);
      }
    }

    out[0][jloc] = cloud_top_pressure[jloc];

    // -----------------------------------------------------------------------------
    // In this section, we calculate the weighted standard deviation about the cloud top pressure.
    double sum_weights = 0.0;
    double weighted_sum_squares = 0.0;
    for (int ilev = nlevs-1; ilev >= tropopauseLevel[jloc]; --ilev) {
      sum_weights += stability_weights[ilev];
      weighted_sum_squares += stability_weights[ilev] *
                                                std::pow(gv_p[ilev] - cloud_top_pressure[jloc], 2);
    }
    standard_deviation[jloc] = std::sqrt(weighted_sum_squares / sum_weights);

    // -----------------------------------------------------------------------------
    // In this section, we find the overcast brightness temperature at the true cloud top pressure
    // by using interpolation.

    // Find the model levels closest to the true cloud top pressure.
    int cloud_model_level = nlevs-1;
    while (gv_p[cloud_model_level] > cloud_top_pressure[jloc]) {
      cloud_model_level--;
    }

    // Interpolate the overcast brightness temperature at the true cloud top pressure.
    // If cloud_model_level == nlevs-1 then the cloud top pressure falls at the lowest model level.
    // It is not possible for the cloud top pressure to fall below the lowest model level because
    // of the way it is assigned above.
    // If cloud_model_level is any other value, a linear interpolation is performed by considering
    // the two model levels adjacent to the cloud top pressure.
    if (cloud_model_level == nlevs-1) {
      bt_cloud_top[jloc] = overcast_bt[cloud_model_level][jloc];
    } else {
      const double bt_diff = overcast_bt[cloud_model_level][jloc] -
                                                    overcast_bt[cloud_model_level + 1][jloc];
      const double log_p_diff = gv_log_p[cloud_model_level] - gv_log_p[cloud_model_level + 1];
      const double bt_gradient = bt_diff / log_p_diff;
      bt_cloud_top[jloc] = overcast_bt[cloud_model_level + 1][jloc]
                      + bt_gradient * (std::log(cloud_top_pressure[jloc]) -
                                          gv_log_p[cloud_model_level + 1]);
    }
  }

  // Output into ObsSpace.
  obsdb.put_db("MetaData", "overcastBriTempAtCloudTop", bt_cloud_top);
  obsdb.put_db("MetaData", "standardDeviationAtCloudTop", standard_deviation);
  obsdb.put_db("DiagnosticFlags/NoStableLayers", "brightnessTemperature", nostablelayers_flag);
  oops::Log::trace() << "StableLayersCloudTopPressure compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

// Function to calculate saturated adiabatic lapse rate.
double StableLayersCloudTopPressure::computeSALR(const double p, const double t) const {
    // Use the Clausius-Clapeyron equation to calculate the saturation vapour pressure.
    const double es = Constants::es_w_0 * exp((Constants::L_c / Constants::rv) *
                                  ((1.0 / Constants::t0c) - (1.0 / t)));

    // Calculate the saturation mixing ratio.
    const double rs = (Constants::epsilon * es) / (p - es);

    // Calculate the saturated adiabatic lapse rate.
    const double numerator = (Constants::rd * t) / (Constants::cp) +
                                      (Constants::L_c * rs) / (Constants::cp);
    const double denominator = p * (1.0 + (std::pow(Constants::L_c, 2) * rs * Constants::epsilon) /
                        (Constants::cp * Constants::rd * std::pow(t, 2)));
    return numerator / denominator;
}

// -----------------------------------------------------------------------------

const ufo::Variables & StableLayersCloudTopPressure::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}   // namespace ufo
