/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Author: Rong Kong (NCAR/MMM) - C++ rewrite of the PPRO operator from Fortran.
 */

#ifndef UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_
#define UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/radarpolarimetric/ppro/Zhang21Forward.h"
#include "ufo/operators/radarshared/MicrophysicsOptions.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

  // PPRO-specific operator type enumeration
  enum class PolarimetricOperatorOption {
    ZHANG21, TCWA2
  };

  struct PolarimetricOperatorOptionParameterTraitsHelper {
    typedef PolarimetricOperatorOption EnumType;
    static constexpr char enumTypeName[] = "PolarimetricOperatorOption";
    static constexpr util::NamedEnumerator<PolarimetricOperatorOption> namedValues[] = {
      { PolarimetricOperatorOption::ZHANG21, "Zhang21" },
      { PolarimetricOperatorOption::TCWA2, "TCWA2" }
    };
  };
}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::PolarimetricOperatorOption> :
    public EnumParameterTraits<ufo::PolarimetricOperatorOptionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------
class ObsPPROParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsPPROParameters, ObsOperatorParametersBase )

 public:
  /// Vertical Coordinate
  oops::Parameter<std::string> VertCoord
    {"VertCoord",
     "vertical coordinate used by the model",
     // this should be consistent with var_prs defined in ufo_vars_mod
     "height_above_mean_sea_level",
     this};

  oops::Parameter<MicrophysicsOption> microOption
    {"microphysics option",
     "microphysics option by name (e.g. Thompson)",
     MicrophysicsOption::THOMPSON,
     this};

  oops::Parameter<PolarimetricOperatorOption> pproOption
    {"polarimetric operator",
     "polarimetric operator by name (Zhang21 or TCWA2)",
     PolarimetricOperatorOption::ZHANG21,
     this};

  oops::Parameter<double> coeff_melt{
    "tuning coefficient for melting",                // YAML name
    "Coefficient scaling the geometric mean sqrt(qr*qx)/sqrt(nr*qx) "
    "in Liu et al. (2024) melting model",            // description
    0.3,                                             // default
    this
  };

  oops::Parameter<std::string> melting_scheme{
    "melting scheme",                                // YAML name
    "Melting scheme option: 'liu24' (default) or 'jung08'.",
    "liu24",                                         // default: liu24
    this
  };

  // Power-law exponents in the liu24 melting model:
  //   qmx  = qr^(melting exponent rain) * qx^(melting exponent ice)
  //   ratx = qr^(fraction exponent rain)
  //          / (qr^(fraction exponent rain) + qx^(fraction exponent ice))
  // Defaults reproduce the geometric mean sqrt(qr*qx) and standard qr/(qr+qx).
  oops::Parameter<double> melting_exp_rain{
    "melting exponent rain",
    "Exponent on qr / ntr in the liu24 melting-water and number-concentration "
    "geometric mean (default 0.5).",
    0.5, this
  };
  oops::Parameter<double> melting_exp_ice{
    "melting exponent ice",
    "Exponent on qx / ntx in the liu24 melting-water and number-concentration "
    "geometric mean (default 0.5).",
    0.5, this
  };
  oops::Parameter<double> fraction_exp_rain{
    "fraction exponent rain",
    "Exponent on qr in the liu24 melting liquid-fraction ratio (default 1.0).",
    1.0, this
  };
  oops::Parameter<double> fraction_exp_ice{
    "fraction exponent ice",
    "Exponent on qx in the liu24 melting liquid-fraction ratio (default 1.0).",
    1.0, this
  };

  oops::Parameter<bool> enable_melting_transition{
    "enable melting transition",                     // YAML name
    "Enable smooth transition function for melting based on qx/qr ratio balance",  // description
    false,                                           // default: disabled
    this
  };

  oops::Parameter<double> snow_ratio_low{
    "snow ratio low",
    "Lower bound of transition for snow (factor=0 below this)",
    0.01,
    this
  };

  oops::Parameter<double> snow_ratio_high{
    "snow ratio high",
    "Upper bound of transition for snow (factor=1 above this)",
    0.05,
    this
  };

  oops::Parameter<double> graupel_ratio_low{
    "graupel ratio low",
    "Lower bound of transition for graupel (factor=0 below this)",
    0.01,
    this
  };

  oops::Parameter<double> graupel_ratio_high{
    "graupel ratio high",
    "Upper bound of transition for graupel (factor=1 above this)",
    0.05,
    this
  };

  oops::Parameter<double> hail_ratio_low{
    "hail ratio low",
    "Lower bound of transition for hail (factor=0 below this)",
    0.01,
    this
  };

  oops::Parameter<double> hail_ratio_high{
    "hail ratio high",
    "Upper bound of transition for hail (factor=1 above this)",
    0.05,
    this
  };

  // Melting water content limit: qmxr <= melting_water_fraction * qr
  // Default 1.0 = ensure qmxr <= qr (physical constraint);
  // < 1.0 = tighter limit (e.g. 0.3 = max 30%)
  oops::Parameter<double> melting_water_fraction{
    "melting water fraction",
    "Max fraction of qr allowed for melting water content (qmsr/qmgr/qmhr). "
    "1.0 (default) = qmxr <= qr (basic physical constraint); "
    "0.3 = qmxr <= 0.3*qr (tighter limit, 30% max).",
    1.0,
    this
  };

  // dmmax configuration: maximum mean diameter limits for hydrometeor species [mm]
  // Each species can be individually configured. If not specified, uses default.
  oops::Parameter<double> dmmax_rain{
    "dmmax rain",
    "Max mean diameter for rain [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_snow{
    "dmmax pure snow",
    "Max mean diameter for pure snow [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_snow{
    "dmmax melting snow",
    "Max mean diameter for melting snow [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_graupel{
    "dmmax pure graupel",
    "Max mean diameter for pure graupel [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_graupel{
    "dmmax melting graupel",
    "Max mean diameter for melting graupel [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_hail{
    "dmmax pure hail",
    "Max mean diameter for pure hail [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_hail{
    "dmmax melting hail",
    "Max mean diameter for melting hail [mm]",
    5.0,
    this
  };

  oops::Parameter<bool> treat_hail_as_graupel{
    "treat hail as graupel",
    "When true, use graupel coefficients and formula for hail calculations "
    "(avoids C-band specific hail formula)",
    false,
    this
  };

  // Dm-based melting limit parameters
  oops::Parameter<bool> enable_dm_melting_limit{
    "enable dm melting limit",
    "Enable temperature-dependent Dm-based melting inhibition. When enabled, "
    "large pure ice particles (Dm > threshold) are less likely to melt. "
    "Threshold depends on temperature: 1.0 mm at T < -5°C, 2.0 mm at T > 0°C, "
    "linear interpolation in between.",
    false,  // Default: disabled
    this
  };

  oops::Parameter<double> dm_melting_transition_width{
    "dm melting transition width",
    "Width of smooth transition region for Dm-based melting inhibition [mm]. "
    "Controls how smoothly the melting factor transitions from 1.0 to 0.0 "
    "around the Dm threshold.",
    0.5,  // Default: 0.5 mm
    this
  };

  // ====================================================================
  // Dm tuning coefficients for each hydrometeor type.
  // Multiplicative scaling factors applied to mean diameter (Dm) after computation:
  //   Dm_new = factor * Dm_old,  Z_new = Z_old * factor^3  (since Z ∝ Dm^3)
  // Default 1.0 (no scaling). Values <1 reduce Dm, >1 increase it.
  // Each coefficient affects all downstream radar variables (ZH, ZDR, KDP)
  // for the corresponding species.
  // ====================================================================
  oops::Parameter<double> tuning_dm_rain{
    "tuning dm rain",
    "Multiplicative scaling factor for rain mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_melting_snow{
    "tuning dm melting snow",
    "Multiplicative scaling factor for melting snow mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_melting_graupel{
    "tuning dm melting graupel",
    "Multiplicative scaling factor for melting graupel mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_melting_hail{
    "tuning dm melting hail",
    "Multiplicative scaling factor for melting hail mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_pure_snow{
    "tuning dm pure snow",
    "Multiplicative scaling factor for pure snow mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_pure_graupel{
    "tuning dm pure graupel",
    "Multiplicative scaling factor for pure graupel mean diameter (Dm).",
    1.0, this
  };
  oops::Parameter<double> tuning_dm_pure_hail{
    "tuning dm pure hail",
    "Multiplicative scaling factor for pure hail mean diameter (Dm).",
    1.0, this
  };

  // ====================================================================
  // Melting fraction tuning coefficients.
  // Multiplicative scaling factors applied to the melting fraction
  // (ratio = qr/(qr+qx)) after it is computed by the melting scheme:
  //   ratio_new = min(factor * ratio_old, 1.0)
  // Default 1.0 (no scaling). <1 makes particles more ice-like,
  // >1 makes them more rain-like. Clamped to [0, 1].
  // Affects: melting density, scattering coefficients, mass partitioning.
  // ====================================================================
  oops::Parameter<double> tuning_melt_frac_snow{
    "tuning melt frac snow",
    "Multiplicative scaling factor for melting snow liquid fraction (rats = qr/(qr+qs)).",
    1.0, this
  };
  oops::Parameter<double> tuning_melt_frac_graupel{
    "tuning melt frac graupel",
    "Multiplicative scaling factor for melting graupel liquid fraction (ratg = qr/(qr+qg)).",
    1.0, this
  };
  oops::Parameter<double> tuning_melt_frac_hail{
    "tuning melt frac hail",
    "Multiplicative scaling factor for melting hail liquid fraction (rath = qr/(qr+qh)).",
    1.0, this
  };

  // Skip small qx: when enabled, skips all radar variable calculations for species
  // with negligible mixing ratios, setting ZH/KDP/dm/w/z=0, ZDR=1, PhV=1.
  oops::Parameter<bool> skip_small_qx{
    "skip small qx",
    "Skip radar variable calculations for hydrometeor species with very small mixing ratios "
    "(qpr < 1e-3, qms < 0.001, qmg < 0.01, qmh < 0.01, qps < 1e-4, qpg < 0.001, qph < 0.001 g/kg). "
    "Sets ZH/KDP/dm=0, ZDR=1.0 (0 dB), PhV=1.0 for skipped species.",
    true,  // Default: enabled
    this
  };

  /*
   * The following list of hydrometeor species mixing ratios and number concentrations
   * need to be consistent with list of variables in ufo_variables_mod, which has some
   * confusing names that need to be resolved in the future - all to match CCPP
   * convention.  In the meantime, FV3 and MPAS use different variables.  The code
   * will default to FV3 names (rain_water, snow_water, and graupel), but all input
   * species can be overridden with these optional YAML parameters, which for MPAS
   * is mass_content_of_rain_in_atmosphere_layer for example.
   */

  oops::Parameter<std::string> var_rain_mixing_ratio
    {"var_rain_mixing_ratio",
     "Name of model rain mixing ratio variable",
     "rain_water",  // this should be consistent with var_qr
     this};

  oops::Parameter<std::string> var_snow_mixing_ratio
    {"var_snow_mixing_ratio",
     "Name of model snow mixing ratio variable",
     "snow_water",  // this should be consistent with var_qs
     this};

  oops::Parameter<std::string> var_graupel_mixing_ratio
    {"var_graupel_mixing_ratio",
     "Name of model graupel mixing ratio variable",
     "graupel",  // this should be consistent with var_qg
     this};

  oops::Parameter<std::string> var_hail_mixing_ratio
    {"var_hail_mixing_ratio",
     "Name of model hail mixing ratio variable",
     "hail",  // this should be consistent with var_qh
     this};

  oops::Parameter<std::string> var_rain_number_concentration
    {"var_rain_number_concentration",
     "Name of model rain number concentration variable",
     "rain_number_concentration",  // this should be consistent with var_nr
     this};

  oops::Parameter<std::string> var_snow_number_concentration
    {"var_snow_number_concentration",
     "Name of model snow number concentration variable",
     "snow_number_concentration",  // this should be consistent with var_ns
     this};

  oops::Parameter<std::string> var_graupel_number_concentration
    {"var_graupel_number_concentration",
     "Name of model graupel number concentration variable",
     "graupel_number_concentration",  // this should be consistent with var_ng
     this};

  oops::Parameter<std::string> var_hail_number_concentration
    {"var_hail_number_concentration",
     "Name of model hail number concentration variable",
     "hail_number_concentration",  // this should be consistent with var_nh
     this};

  oops::Parameter<std::string> var_graupel_vol_mixing_ratio
    {"var_graupel_vol_mixing_ratio",
     "Name of model graupel volume mixing ratio variable",
     "volume_mixing_ratio_of_graupel_in_air",  // this should be consistent with var_qvg
     this};

  oops::Parameter<std::string> var_hail_vol_mixing_ratio
    {"var_hail_vol_mixing_ratio",
     "Name of model hail volume mixing ratio variable",
     "volume_mixing_ratio_of_hail_in_air",  // this should be consistent with var_qvh
     this};

  // TCWA2 operator specific variables (NOT needed for Zhang21)
  // These are only used when polarimetric_operator = TCWA2
  oops::Parameter<std::string> var_cloud_mixing_ratio
    {"var_cloud_mixing_ratio",
     "Name of model cloud water mixing ratio variable "
     "(TCWA2 operator only, not needed for Zhang21)",
     "cloud_water",  // this should be consistent with var_qc
     this};

  oops::Parameter<std::string> var_ice_mixing_ratio
    {"var_ice_mixing_ratio",
     "Name of model ice mixing ratio variable (TCWA2 operator only, not needed for Zhang21)",
     "ice_water",  // this should be consistent with var_qi
     this};

  oops::Parameter<std::string> var_ice_number_concentration
    {"var_ice_number_concentration",
     "Name of model ice number concentration variable "
     "(TCWA2 operator only, not needed for Zhang21)",
     "ice_number_concentration",  // this should be consistent with var_ni
     this};

  oops::Parameter<std::string> var_snow_melted_fraction
    {"var_snow_melted_fraction",
     "Name of model snow melted fraction variable (TCWA2 operator only, not needed for Zhang21)",
     "snow_melted_fraction",  // this should be consistent with var_smlf
     this};

  oops::Parameter<std::string> var_graupel_melted_fraction
    {"var_graupel_melted_fraction",
     "Name of model graupel melted fraction variable (TCWA2 operator only, not needed for Zhang21)",
     "graupel_melted_fraction",  // this should be consistent with var_gmlf
     this};

  // Dm regularization for small ntx/qx
  oops::Parameter<bool> enable_dm_regularization{
    "enable dm regularization",
    "Enable Dm regularization for small ntx/qx in dm_z_2moment. When enabled, "
    "sets dm = 0 if both qx < 1e-3 g/kg AND ntx < 1000 /m^3 (negligible "
    "particle population). Otherwise uses standard formula (caller ensures "
    "ntx >= 10 for zero protection, dmmax caps upper bound). "
    "Default: true (enabled).",
    true,   // Default: enabled
    this
  };
};

/// RadarReflectivity observation operator class
class ObsPPRO : public ObsOperatorBase,
                   private util::ObjectCounter<ObsPPRO> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsPPROParameters Parameters_;
  static const std::string classname() {return "ufo::ObsPPRO";}

  ObsPPRO(const ioda::ObsSpace &, const ObsPPROParameters &);
  virtual ~ObsPPRO();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                  const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;

  /// Build the Zhang21 tuning parameters from the YAML-derived operator params.
  static ppro::PPROTuning buildTuning(const ObsPPROParameters &params);

  const ioda::ObsSpace& odb_;
  oops::Variables varin_;

  /// Microphysics scheme selecting which model variables are required.
  MicrophysicsOption microOption_;
  /// True if the TCWA2 forward operator is selected (else Zhang21).
  bool useTcwa2_;
  /// Melting coefficient passed to Zhang21 (NSSL/2-moment-with-hail only).
  double coeffMelt_;
  /// Name of the model vertical-coordinate GeoVaL used for interpolation.
  std::string vCoord_;
  /// Ordered list of model variable names required from GeoVaLs (excl. vCoord_).
  std::vector<std::string> geovars_;
  /// ObsVector column index for each operator output (zh, zdr, kdp, phv);
  /// -1 if the corresponding variable is not simulated.
  std::array<int, 4> outColumn_;

  /// Zhang21 forward operator (holds precomputed coefficients + tuning).
  ppro::Zhang21Forward zhang21_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_
