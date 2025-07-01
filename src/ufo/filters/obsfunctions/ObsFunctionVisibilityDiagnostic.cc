/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionVisibilityDiagnostic.h"

#include <algorithm>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

static ObsFunctionMaker<VisibilityDiagnostic> makerVisibilityDiagnostic_(
    "VisibilityDiagnostic");

VisibilityDiagnostic::VisibilityDiagnostic(
    const eckit::LocalConfiguration& conf) {
  oops::Log::trace() << "VisibilityDiagnostic constructor" << std::endl;
  // Check options
  options_.validateAndDeserialize(conf);
}

VisibilityDiagnostic::~VisibilityDiagnostic() {
  oops::Log::trace() << "VisibilityDiagnostic destructor" << std::endl;
}

void VisibilityDiagnostic::compute(const ObsFilterData& in,
                                   ioda::ObsDataVector<float>& out) const {
  oops::Log::trace() << "VisibilityDiagnostic compute start" << std::endl;
  if (out.nvars() != 1) {
    throw eckit::BadValue("VisibilityDiagnostic: one output variable expected",
                          Here());
  }
  const float missing = util::missingValue<float>();
  const size_t nlocs = in.nlocs();

  // Input Variables (see header file for details)
  std::vector<float> rh(nlocs), p(nlocs), t(nlocs);
  in.get(options_.rh.value(), rh);
  in.get(options_.p.value(), p);
  in.get(options_.t.value(), t);

  // Controlling parameters (see header file for details)
  const float rh_crit = options_.rh_crit.value();
  const float ccp1 = options_.ccp1.value();
  const float ccp2 = options_.ccp2.value();
  const float ccp3 = options_.ccp3.value();
  const float rh_tot_p1 = options_.rh_tot_p1.value();
  const float rh_tot_p2 = options_.rh_tot_p2.value();
  const float rh_tot_p3 = options_.rh_tot_p3.value();
  const float m = options_.m.value();
  const float rho_aerosol = options_.rho_aerosol.value();
  const float r_0 = options_.r_0.value();
  const float N_0 = options_.N_0.value();
  const float rho_air = options_.rho_air.value();
  const float rho_wat = options_.rho_wat.value();
  const float size_distribution_exponent =
      options_.size_distribution_exponent.value();
  const float A = options_.A.value();
  const float B = options_.B.value();
  const float liminal_contrast = options_.liminal_contrast.value();
  const float extinction_efficiency = options_.extinction_efficiency.value();
  const float particle_size_weighting =
      options_.particle_size_weighting.value();
  const float visibility_limit = options_.visibility_limit.value();
  const float q_tot_min = options_.q_tot_min.value();
  const float rh_tot_min = options_.rh_tot_min.value();
  const float rh_tot_max = options_.rh_tot_max.value();
  const float g_min = options_.g_min.value();
  const float g_max = options_.g_max.value();
  const int max_iterations = options_.max_iterations.value();
  const float max_inflation_factor = options_.max_inflation_factor.value();
  const float min_inflation_factor = options_.min_inflation_factor.value();
  const float delta_g = options_.delta_g.value();
  const float alpha = options_.alpha.value();

  // Derived values
  constexpr auto clark2008 = formulas::Formulation::Clark2008;
  const float m_0 = formulas::dryAerosolMassMixingRatio(rho_aerosol, rho_air,
                                                        N_0, r_0, clark2008);
  const float r_md = formulas::dryMeanParticleVolumeRadius(
      r_0, m, m_0, size_distribution_exponent, clark2008);
  const float N = formulas::aerosolNumberDensity(
      N_0, m, m_0, size_distribution_exponent, clark2008);
  const float r_act =
      formulas::approximateActivationDropletRadius(A, B, r_md, clark2008);
  const float beta_0 = formulas::approximateAerosolExtinctionCoefficientFactor(
      extinction_efficiency, particle_size_weighting, clark2008);

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (rh[iloc] == missing || p[iloc] == missing || t[iloc] == missing) {
      out[0][iloc] = missing;
      continue;
    }
    // Cloud cover is a proxy for the presence of cloud or fog.
    const float cc = formulas::cloudCoverFraction(
        rh[iloc], rh_crit, ccp1, ccp2, ccp3,
        formulas::Formulation::UKMOUnififiedModelOPS);
    const float rh_tot = formulas::totalWaterRelativeHumidity(
        cc, rh_crit, rh_tot_p1, rh_tot_p2, rh_tot_p3,
        formulas::Formulation::UKMOUnififiedModelOPS);
    const float q_sat = saturationSpecificHumidity(t[iloc], p[iloc]);
    const float q_tot =
        totalWaterSpecificHumidity(rh[iloc], rh_crit, rh_tot, q_sat);
    const float r_m = solveForDropletRadius(
        r_act, r_md, q_tot, q_tot_min, rh_tot, rh_tot_min, rh_tot_max, q_sat, A,
        B, N, rho_wat, g_min, g_max, max_iterations, max_inflation_factor,
        min_inflation_factor, delta_g, alpha);
    const float beta =
        formulas::aerosolExtinctionCoefficient(beta_0, N, r_m, clark2008);
    const float visibility =
        formulas::aerosolVisibility(liminal_contrast, beta, clark2008);
    // Limit visibility to avoid unrealistically high values
    out[0][iloc] = 1.0f / (1.0f / visibility + 1.0f / visibility_limit);
  }  // nlocs
  oops::Log::trace() << "VisibilityDiagnostic compute complete" << std::endl;
}

float VisibilityDiagnostic::saturationSpecificHumidity(const float t,
                                                       const float p) const {
  const float esat_pure_water = formulas::SatVaporPres_fromTemp(
      t, formulas::Formulation::GoffGratchLandoltBornsteinWater);
  const float esat_air = formulas::SatVaporPres_correction(
      esat_pure_water, t, p, formulas::Formulation::Gill);
  return formulas::Qsat_From_Psat(esat_air, p, formulas::Formulation::GillUKMO);
}

float VisibilityDiagnostic::totalWaterSpecificHumidity(
    const float rh, const float rh_crit, const float rh_tot,
    const float q_sat) const {
  if (rh <= rh_crit) {
    // No cloud/fog (q_tot = q where q is specific humidity)
    return rh * q_sat;
  } else {
    // Includes cloud/fog, hence using rh_tot
    return rh_tot * q_sat;
  }
}

float VisibilityDiagnostic::solveForDropletRadius(
    const float r_act, const float r_md, float q_tot, const float q_tot_min,
    float rh_tot, const float rh_tot_min, const float rh_tot_max,
    const float q_sat, const float A, const float B, const float N,
    const float rho_wat, const float g_min, const float g_max,
    const int max_iterations, const float max_inflation_factor,
    const float min_inflation_factor, const float delta_g,
    const float alpha) const {
  // The droplet radius is solved for in terms of radii normalized by the
  // dry mean particle volume radius, so the activation droplet radius also
  // needs normalizing.
  const float g_act = r_act / r_md;
  // Limit total water:
  q_tot = std::max(q_tot, q_tot_min);
  rh_tot = std::max(q_tot / q_sat, rh_tot_min);
  rh_tot = std::min(rh_tot, rh_tot_max);
  // Calculate the first guess g_0 by taking Eq. 10 in Clark et al. (see header
  // for full reference) and assuming r_m > r_md (i.e. g > 1), then rearrange
  // for g = r_m / r_md
  const float g_0 = std::pow(1.0f - B / std::log(rh_tot), 1.0f / 3.0f);
  // Current guess of g
  float g_n = g_0;
  // Next guess of g
  float g_n_p_1 = g_n;
  for (int n = 0; n < max_iterations; ++n) {
    // Exponential smoothing
    g_n = alpha * g_n_p_1 + (1 - alpha) * g_n;
    g_n_p_1 =
        newtonRaphsonStep(g_n, g_act, r_md, q_tot, q_sat, A, B, N, rho_wat);
    g_n_p_1 = std::min(g_max, g_n_p_1);
    g_n_p_1 = std::max(g_min, g_n_p_1);
    g_n_p_1 = std::min(g_n_p_1, g_n * max_inflation_factor);
    g_n_p_1 = std::max(g_n_p_1, g_n * min_inflation_factor);
    if (std::abs(g_n_p_1 - g_n) < delta_g) {
      break;
    }
  }
  // Get the droplet radius r_m (m)
  return r_md * g_n_p_1;
}

float VisibilityDiagnostic::newtonRaphsonStep(
    const float g_n, const float g_act, const float r_md, const float q_tot,
    const float q_sat, const float A, const float B, const float N,
    const float rho_wat) const {
  // Perform Newton-Raphson step, i.e. aiming to have f = 0.0
  float f;
  float f_prime;
  constexpr auto clark2008 = formulas::Formulation::Clark2008;
  const float q_l =
      formulas::specificCloudWaterContent(g_n, r_md, N, rho_wat, clark2008);
  const float q_l_prime = formulas::specificCloudWaterContentDerivative(
      g_n, r_md, N, rho_wat, clark2008);
  if (g_n < g_act) {
    // Below peak of Koehler curve
    f = q_tot -
        q_sat *
            formulas::equilibriumRelativeHumidity(g_n, r_md, A, B, clark2008) -
        q_l;
    f_prime = -q_sat * formulas::equilibriumRelativeHumidityDerivative(
                           g_n, r_md, A, B, clark2008) -
              q_l_prime;
  } else {
    // At or beyond peak of Koehler curve - force f back towards g_act (where
    // equilibriumRelativeHumidityDerivative = 0.0)
    f = q_tot -
        q_sat * formulas::equilibriumRelativeHumidity(g_act, r_md, A, B,
                                                      clark2008) -
        q_l;
    f_prime = -q_l_prime;
  }
  return g_n - f / f_prime;
}

const ufo::Variables& VisibilityDiagnostic::requiredVariables() const {
  return invars_;
}

}  // namespace ufo
