/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVISIBILITYDIAGNOSTIC_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVISIBILITYDIAGNOSTIC_H_

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/variabletransforms/Cal_Humidity.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsFunctionVisibilityDiagnostic ObsFunction
class VisibilityDiagnosticParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VisibilityDiagnosticParameters, Parameters)

 public:
  /// Relative humidity units for the input relative humidity variable and other
  /// relative humidity parameters, either percentage or fraction.
  oops::RequiredParameter<RelativeHumidityUnits> rhUnits{
      "observation relative humidity units", this};

  /// The below parameters are required observed or derived variables

  /// Relative humidity (% or fraction) - defined as the ratio of the specific
  /// humidity to the saturation specific humidity
  oops::RequiredParameter<Variable> rh{"relative humidity variable", this};
  /// Pressure (Pa)
  oops::RequiredParameter<Variable> p{"pressure variable", this};
  /// Temperature (K)
  oops::RequiredParameter<Variable> t{"temperature variable", this};

  /// Critical relative humidity (unitless) - that at which liquid water
  /// droplets are considered to form (fraction, not percentage)
  oops::RequiredParameter<float> rh_crit{
      "critical relative humidity",
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};

  /// The below parameterize the `RH_TO_CC` diagnostic taken from the UK Met
  /// Office Unified Model, as used in the Met Office Observation Processing
  /// System (OPS). This calculates a cloud cover fraction, a proxy for the
  /// presence of mist or fog, from relative humidity.

  /// Cloud cover parameter 1 (unitless) - used to estimate liquid water content
  oops::RequiredParameter<float> ccp1{"cloud cover parameter 1", this};
  /// Cloud cover parameter 2 (unitless) - used to estimate liquid water content
  oops::RequiredParameter<float> ccp2{"cloud cover parameter 2", this};
  /// Cloud cover parameter 3 (unitless) - used to estimate liquid water content
  oops::RequiredParameter<float> ccp3{"cloud cover parameter 3", this};

  /// The below parameterize the CC_TO_RHTOT diagnostic taken from the UK Met
  /// Office Unified Model, as used in the Met Office Observation Processing
  /// System (OPS). In that diagnostic, total water relative humidity is
  /// associated with cloud cover via a piecewise function, derived from a
  /// triangular distribution. This is consistent with the Smith (1990) cloud
  /// scheme (https://doi.org/10.1002/qj.49711649210).

  /// Total water relative humidity parameter 1 - The cloud cover fraction at
  /// which the piecewise function transitions.
  oops::RequiredParameter<float> rh_tot_p1{
      "total water relative humidity parameter 1", this};
  /// Total water relative humidity parameter 2 - The scaling to apply to the
  /// cloud cover fraction within the first part of the piecewise function.
  oops::RequiredParameter<float> rh_tot_p2{
      "total water relative humidity parameter 2", this};
  /// Total water relative humidity parameter 3 - The scaling to apply to (1 -
  /// cc), where cc is the cloud cover fraction, within the second part of the
  /// piecewise function.
  oops::RequiredParameter<float> rh_tot_p3{
      "total water relative humidity parameter 3", this};

  /// The below parameters correspond to variables found in Clark, Peter A.,
  /// et al. "Prediction of visibility and aerosol within the operational Met
  // Office Unified Model. I: Model formulation and variational assimilation."
  /// Quarterly Journal of the Royal Meteorological Society 134.636 (2008):
  /// 1801-1816. DOI 10.1002/qj.318. The symbol names which correspond to these
  /// parameters in the paper are in backticks as latex variables.

  /// Aerosol mass concentration `m` in Eq. (1) (kg kg^-1) - also known as the
  /// dry aerosol mass mixing ratio
  oops::RequiredParameter<float> m{"aerosol mass concentration",
                                   this,
                                   {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Aerosol density `\rho` in Eq. 1 (kg m^-3)
  oops::Parameter<float> rho_aerosol{
      "aerosol density",
      1700.0f,
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Standard dry aerosol radius `r_0` in Eq. (2) (m)
  oops::Parameter<float> r_0{"standard aerosol radius",
                             0.16E-6f,
                             this,
                             {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Standard aerosol number density `N_0` in Eq. (3) (m^-3)
  oops::Parameter<float> N_0{"standard aerosol number density",
                             5.E8f,
                             this,
                             {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Air density `\rho_\text{a}` in Eq. (1) (kg m^-3)
  oops::Parameter<float> rho_air{
      "air density", 1.0f, this, {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Water density `\rho_\text{wat}` in Eq. (13) (kg m^-3)
  oops::Parameter<float> rho_wat{"water density",
                                 1000.0f,
                                 this,
                                 {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Aerosol size distribution power law exponent `p` in Eq. (2)
  oops::Parameter<float> size_distribution_exponent{
      "aerosol size distribution power law exponent", 1.0f / 6.0f, this};
  /// Koehler curve constant `A` in Eq. 10, related to the surface tension of
  /// water (m)
  oops::Parameter<float> A{"Koehler curve constant A",
                           1.2E-9f,
                           this,
                           {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Koehler curve constant `B` in Eq. 10, the "activation parameter"
  /// (unitless)
  oops::Parameter<float> B{"Koehler curve constant B", 0.5f, this};
  /// Liminal contrast factor `\epsilon` in Eq. 7 (unitless) for calculation of
  /// extinction coefficient
  oops::Parameter<float> liminal_contrast{
      "liminal contrast",
      0.02f,
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Extinction efficiency factor `Q` in Eq. 16 (unitless) for calculation of
  /// extinction coefficient
  oops::Parameter<float> extinction_efficiency{"extinction efficiency factor Q",
                                               2.0f, this};
  /// Particle size weighting factor `\eta` in Eq. 19 (unitless) for calculation
  /// of extinction coefficient
  oops::Parameter<float> particle_size_weighting{
      "particle size weighting factor eta", 0.75f, this};
  /// Clean air visibility (m) - used in place of the extinction coefficient of
  /// air to limit unrealistically high visibilities (see text after Eq. 8)
  oops::Parameter<float> visibility_limit{
      "visibility limit",
      100000.0f,
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};

  /// The below parameters are specific to the Newton-Raphson method used to
  /// estimate the particle radius `r_m` (m) given the atmospheric conditions.

  /// Minimum allowed total-water specific humidity `q_T^\text{min}` (kg kg^-1)
  oops::Parameter<float> q_tot_min{"minimum total water specific humidity",
                                   0.001f,
                                   this,
                                   {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Minimum allowed total water relative humidity for calculating
  /// the first guess of the particle radius `r` (fraction or percentage).
  /// Default is 0.01 if `observation relative humidity units` is fraction or
  /// 1.0 % if `observation relative humidity units` is percentage.
  oops::OptionalParameter<float> rh_tot_min{
      "minimum relative humidity for first guess",
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Maximum allowed relative humidity for calculating the first
  /// guess of the particle radius `r` (fraction or percentage). Default is
  /// 0.999 if `observation relative humidity units` is fraction or 99.9 % if
  /// `observation relative humidity units` is percentage.
  oops::OptionalParameter<float> rh_tot_max{
      "maximum relative humidity for first guess",
      this,
      {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Minimum allowed wet particle (droplet) radius, normalized by the mean dry
  /// radius `r_m/r_{md}` (m m^-1) - the expression for `r_{md}` is Eq. 2.
  oops::Parameter<float> g_min{"minimum allowed normalized droplet radius",
                               2.0f,
                               this,
                               {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Maximum allowed wet particle (droplet) radius, normalized by the mean dry
  /// radius `r_m/r_{md}` (m m^-1) - the expression for `r_{md}` is Eq. 2.
  oops::Parameter<float> g_max{"maximum allowed normalized droplet radius",
                               10000.0f,
                               this,
                               {oops::exclusiveMinConstraint<float>(0.0f)}};
  /// Maximum allowed number of iterations for the Newton-Raphson method loop
  oops::Parameter<int> max_iterations{
      "maximum number of iterations", 10, this, {oops::minConstraint<int>(0)}};
  /// Maximum allowed droplet radius inflation factor for the Newton-Raphson
  /// loop (i.e. max(r_m^{n+1}/r_m^n) for each iteration, used as a ceiling
  /// (m m^-1)
  oops::Parameter<float> max_inflation_factor{"maximum inflation factor", 2.0f,
                                              this};
  /// Minimum allowed droplet radius inflation factor for the Newton-Raphson
  /// loop (i.e. min(r_m^{n+1}/r_m^n} for each iteration, used as a floor)
  /// (m m^-1)
  oops::Parameter<float> min_inflation_factor{"minimum inflation factor", 0.5f,
                                              this};
  /// Stopping value for the Newton-Raphson loop, the maximum allowed
  /// difference between the normalized droplet radius `r_m/r_{md}` at iteration
  /// `n` and `n+1` (m m^-1)
  oops::Parameter<float> delta_g{"loop stopping value", 0.01f, this};
  /// Exponential smoothing factor to apply to previous and current droplet
  /// radius at the start of the Newton-Raphson loop before continuing: if 1.0,
  /// then the previous droplet radius becomes the new droplet radius and no
  /// smoothing is applied. If 0.0, the smoothing is infinite and the loop will
  /// always keep the first guess of the droplet radius.
  oops::Parameter<float> alpha{
      "droplet radius exponential smoothing factor",
      0.9f,
      this,
      {oops::minConstraint<float>(0.0f), oops::maxConstraint<float>(1.0f)}};
};

// -----------------------------------------------------------------------------

/// \brief Outputs a visibility value (m) given the input observed variables
/// and a set of parameters.
class VisibilityDiagnostic : public ObsFunctionBase<float> {
 public:
  /// @brief Constructs a VisibilityDiagnostic object.
  ///
  /// @param conf The configuration for the visibility diagnostic.
  explicit VisibilityDiagnostic(const eckit::LocalConfiguration &conf);

  /// @brief Destructor for the VisibilityDiagnostic class.
  ~VisibilityDiagnostic();

  /// @brief Computes the visibility diagnostic.
  ///
  /// This method computes the visibility diagnostic based on input data
  /// (typically observations or geovals) in observation space and stores the
  /// result in the output data vector.
  ///
  /// This follows the method outlined in Clark, Peter A., et al. "Prediction of
  /// visibility and aerosol within the operational Met Office Unified Model.
  /// I: Model formulation and variational assimilation."Quarterly Journal of
  /// the Royal Meteorological Society 134.636 (2008): 1801-1816.
  ///  DOI 10.1002/qj.318.
  ///
  /// The visibility is assumed to only be limited by (potentially wet) aerosol
  /// particles, subject to some visibility limit.
  ///
  ///
  /// @param in The input data in observation space.
  /// @param out The output data vector to store the computed visibility values.
  void compute(const ObsFilterData &in,
               ioda::ObsDataVector<float> &out) const override;

  /// @brief Returns an empty ufo::Variables object.
  ///
  /// This is a required method needed for compilation but is not used in this
  /// class.
  ///
  /// @return An empty ufo::Variables object.
  const ufo::Variables &requiredVariables() const override;

 private:
  VisibilityDiagnosticParameters options_;
  ufo::Variables invars_;

  /// @brief The relative humidity at saturation, which is either 100.0f for
  /// percentage or 1.0f for fraction, depending on the configuration.
  const float relativeHumidityAtSaturation_;

  /// @brief Calculates the saturation specific humidity (kg kg^-1) as a
  /// function of the temperature and pressure.
  ///
  /// @param t The temperature (K).
  /// @param p The pressure (Pa).
  /// @return The calculated saturation specific humidity (kg kg^-1).
  float saturationSpecificHumidity(const float t, const float p) const;

  /// @brief Calculates the total water specific humidity (kg kg^-1), i.e. the
  /// total mass of water in any phase per unit mass of dry air.
  ///
  /// @note If the relative humidity is less than or equal to the critical
  /// relative humidity, this is equal to the specific humidity, since no liquid
  /// water is considered to be in the atmosphere.
  ///
  /// @param rh The relative humidity (fraction).
  /// @param rh_crit The critical relative humidity: that at which liquid water
  /// droplets are considered to form (fraction).
  /// @param rh_tot The total water relative humidity (unitless).
  /// @param q_sat The saturation specific humidity (kg kg^-1).
  /// @return The calculated total water specific humidity (kg kg^-1).
  float totalWaterSpecificHumidity(const float rh, const float rh_crit,
                                   const float rh_tot, const float q_sat) const;

  /// @brief Solves Eq. 14 of Clark et al. for the droplet radius (m) using the
  /// Newton-Raphson method.
  ///
  /// @param r_act The activation droplet radius (m).
  /// @param r_md The dry aerosol radius (m).
  /// @param q_tot The total water specific humidity (kg kg^-1).
  /// @param q_tot_min The minimum allowed total water specific humidity
  /// (kg kg^-1).
  /// @param rh_tot The total water relative humidity (unitless).
  /// @param rh_tot_min The minimum allowed total water relative humidity
  /// (unitless).
  /// @param rh_tot_max The maximum allowed total water relative humidity
  /// (unitless).
  /// @param q_sat The saturation specific humidity (kg kg^-1).
  /// @param A The Koehler curve constant A (m).
  /// @param B The Koehler curve constant B (unitless).
  /// @param N The number density of the aerosol particles (m^-3).
  /// @param rho_wat The density of water (kg m^-3).
  /// @param g_min The minimum allowed droplet radius, normalized by the dry
  /// aerosol radius (r_md) (m m^-1).
  /// @param g_max The maximum allowed droplet radius, normalized by the dry
  /// aerosol radius (r_md) (m m^-1).
  /// @param max_iterations The maximum allowed number of iterations for the
  /// Newton-Raphson loop.
  /// @param max_inflation_factor The maximum allowed droplet radius inflation
  /// factor for the Newton-Raphson loop (i.e. max(r_m^{n+1}/r_m^n) for each
  /// iteration, used as a ceiling) (m m^-1)
  /// @param min_inflation_factor The minimum allowed droplet radius inflation
  /// factor for the Newton-Raphson loop (i.e. min(r_m^{n+1}/r_m^n} for each
  /// iteration, used as a floor) (m m^-1)
  /// @param delta_g The stopping value for the Newton-Raphson loop, i.e. the
  /// maximum allowed difference between the normalized droplet radius
  /// `r_m/r_{md}` at iteration `n` and `n+1` (m m^-1)
  /// @param alpha The weighting to apply to previous and
  /// current droplet radius at the start of the Newton-Raphson loop before
  /// continuing: if 1.0, then the previous droplet radius becomes the new
  /// droplet radius. If 0.0, then the loop will revert to the first guess of
  /// the droplet radius (unitless).
  /// @return The calculated droplet radius r_m (m).
  float solveForDropletRadius(const float r_act, const float r_md, float q_tot,
                              const float q_tot_min, float rh_tot,
                              const float rh_tot_min, const float rh_tot_max,
                              const float q_sat, const float A, const float B,
                              const float N, const float rho_wat,
                              const float g_min, const float g_max,
                              const int max_iterations,
                              const float max_inflation_factor,
                              const float min_inflation_factor,
                              const float delta_g, const float alpha) const;

  /// @brief Calculates the next Newton-Raphson step for the normalised droplet
  /// radius (m) using the current guess of the normalized droplet radius.
  ///
  /// This method calculates the next Newton-Raphson step for the normalized
  /// droplet radius (m m^-1) using the current guess of the normalized droplet
  /// radius. If the input radius g_n is greater than the activation radius
  /// g_act, the system is forced back to g_act.
  /// g_act is where the peak of the Koehler curve occurs. For g_n < g_act,
  /// wet particles tend to evaporate to reach equilibrium with the atmosphere.
  /// Where g_n > g_act, wet particles tend to grow through condensation to
  /// reach equilibrium. In the latter case, Newton-Raphson iterations would
  /// tend to increase the droplet radius and a solution would not be reached.
  ///
  /// @param g_n The current guess of the normalized droplet radius (m m^-1).
  /// @param g_act The activation droplet radius (m m^-1).
  /// @param r_md The dry aerosol radius (m).
  /// @param q_tot The total water specific humidity (kg kg^-1).
  /// @param q_sat The saturation specific humidity (kg kg^-1).
  /// @param A The Koehler curve constant A (m).
  /// @param B The Koehler curve constant B (unitless).
  /// @param N The number density of the aerosol particles (m^-3).
  /// @param rho_wat The density of water (kg m^-3).
  float newtonRaphsonStep(const float g_n, const float g_act, const float r_md,
                          const float q_tot, const float q_sat, const float A,
                          const float B, const float N,
                          const float rho_wat) const;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVISIBILITYDIAGNOSTIC_H_
