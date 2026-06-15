/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_ZHANG21FORWARD_H_
#define UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_ZHANG21FORWARD_H_

#include <optional>
#include <string>

// ----------------------------------------------------------------------------
// Zhang21 Polarimetric Radar Operator (C++ port of zhang21_forward_mod.f90)
//
// Parameterized polarimetric radar operator following
//   Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for
//   simulation and assimilation of polarimetric radar data with numerical
//   weather predictions. Adv. Atmos. Sci., 38(5), 737-754.
//
// Author: Rong Kong (NCAR/MMM) - C++ rewrite of the operator from Fortran;
//         bug fixes, operator tuning and testing, modularization and
//         optimization.
//
// Original Fortran authors:
//   Zhiquan (Jake) Liu (NCAR/MMM) - Initial implementation, July 2022.
//   Hejun Xie (NCAR/MMM) - Integration into JEDI-UFO framework, August 2024.
//   Tao Sun (NCAR/MMM) - Hail categories and extended microphysics, April 2025.
// ----------------------------------------------------------------------------

namespace ufo {
namespace ppro {

/// Runtime-tunable parameters for the Zhang21 operator (mirrors the module-level
/// `save` variables in zhang21_forward_mod.f90 that are set from YAML).
struct PPROTuning {
  bool   skip_small_qx = true;

  bool   enable_melting_transition = false;
  double snow_ratio_low = 0.01;
  double snow_ratio_high = 0.05;
  double graupel_ratio_low = 0.01;
  double graupel_ratio_high = 0.05;
  double hail_ratio_low = 0.01;
  double hail_ratio_high = 0.05;

  double melting_water_fraction = 1.0;

  bool   enable_dm_melting_limit = false;
  double dm_melting_transition_width = 0.5;

  bool   enable_dm_regularization = true;

  double tuning_dm_rain = 1.0;
  double tuning_dm_melting_snow = 1.0;
  double tuning_dm_melting_graupel = 1.0;
  double tuning_dm_melting_hail = 1.0;
  double tuning_dm_pure_snow = 1.0;
  double tuning_dm_pure_graupel = 1.0;
  double tuning_dm_pure_hail = 1.0;

  double melting_exp_rain = 0.5;
  double melting_exp_ice = 0.5;

  double fraction_exp_rain = 1.0;
  double fraction_exp_ice = 1.0;

  double tuning_melt_frac_snow = 1.0;
  double tuning_melt_frac_graupel = 1.0;
  double tuning_melt_frac_hail = 1.0;

  std::string melting_scheme_option = "liu24";

  double dmmax_rain = 5.0;
  double dmmax_pure_snow = 10.0;
  double dmmax_melting_snow = 5.0;
  double dmmax_pure_graupel = 10.0;
  double dmmax_melting_graupel = 5.0;
  double dmmax_pure_hail = 10.0;
  double dmmax_melting_hail = 5.0;

  bool   treat_hail_as_graupel = false;
};

/// Optional per-observation inputs (mirror the optional arguments of
/// zhang21_compute_point).  A disengaged optional means "not present".
struct Zhang21Optionals {
  std::optional<double> qh;
  std::optional<double> nr;
  std::optional<double> ns;
  std::optional<double> ng;
  std::optional<double> nh;
  std::optional<double> vg;
  std::optional<double> vh;
  std::optional<double> coeff_melt;
  std::optional<double> temperature;  // Celsius
};

class Zhang21Forward {
 public:
  explicit Zhang21Forward(const PPROTuning &tuning);

  /// Compute dual-pol variables (linear units) for a single observation point.
  ///   iband        : 1=S-band, 2=C-band
  ///   density_air  : dry air density [kg/m^3]
  ///   temp_air     : air temperature [K]
  ///   qr,qs,qg     : rain/snow/graupel mixing ratio [g/kg]
  void computePoint(int iband, double density_air, double temp_air,
                    double qr, double qs, double qg,
                    double &zh, double &zdr, double &kdp, double &phv,
                    const Zhang21Optionals &opt) const;

 private:
  // Hydrometeor densities [g/cm^3]
  static constexpr double density_rain    = 1.0;
  static constexpr double density_snow    = 0.1;
  static constexpr double density_graupel = 0.5;
  static constexpr double density_hail    = 0.917;

  static constexpr double qsum_min_threshold = 1.0e-3;  // g/kg

  PPROTuning t_;

  // Precomputed scattering "a" coefficients for dry hydrometeors (ratio=0),
  // i.e. a(i) = iceCoef(coefs, i, 0).  Indexed [0..12].
  double sband_snow_a_[13];
  double sband_graupel_a_[13];
  double sband_hail_a_[13];
  double cband_snow_a_[13];
  double cband_graupel_a_[13];
  double cband_hail_a_[13];

  // Melting-scheme output bundle (mirrors the many output args of the Fortran
  // melting_model_* subroutines).
  struct MeltOut {
    double qms = 0.0, qmg = 0.0, qpr = 0.0, qps = 0.0, qpg = 0.0;
    double rats = 0.0, ratg = 0.0;
    double qmh = 0.0, qph = 0.0, rath = 0.0;
    double ntms = 0.0, ntmg = 0.0, ntpr = 0.0, ntps = 0.0, ntpg = 0.0;
    double ntmh = 0.0, ntph = 0.0;
  };

  // --- Low-level helpers (1:1 with the Fortran subroutines) ------------------
  static double coefA(const double c[4], double gamma);

  static void dualpolOpRain(double w, double dm, const double *rainflat,
                            double &zh, double &zdr, double &kdp, double &phv);

  static void dualpolOpIcephase(double z, double w, double rho, double dm,
                                const double a[13],
                                double &zh, double &zdr, double &kdp, double &phv,
                                int iband, const std::string &ptype);

  static void dualpolOpTotal(double zh_prain, double zh_psnow, double zh_pgraupel,
                             double zdr_prain, double zdr_psnow, double zdr_pgraupel,
                             double kdp_prain, double kdp_psnow, double kdp_pgraupel,
                             double phv_prain, double phv_psnow, double phv_pgraupel,
                             double &zh, double &zdr, double &kdp, double &phv,
                             double zh_msnow, double zh_mgraupel,
                             double zdr_msnow, double zdr_mgraupel,
                             double kdp_msnow, double kdp_mgraupel,
                             double phv_msnow, double phv_mgraupel,
                             bool include_hail,
                             double zh_mhail, double zdr_mhail,
                             double kdp_mhail, double phv_mhail,
                             double zh_phail, double zdr_phail,
                             double kdp_phail, double phv_phail);

  static double weticephaseDensity(int species, double ratio);  // 0=snow,1=graupel,2=hail

  static double waterContent(double air_density, double qx) { return air_density * qx; }

  void dmZ2moment(double air_density, double qx, double rhox, double ntx,
                  double &dm, double &z,
                  std::optional<double> dmmax, std::optional<double> dm_min) const;

  static void n0LambdaWsm6(int precip_type, double air_density, double temp,
                           double qx, double rhox, double &n0, double &lambda);

  static void dmZwsm6(int precip_type, double air_density, double temp,
                      double qx, double rhox, double &dm, double &z,
                      std::optional<double> dmmax, std::optional<double> dm_min);

  static double getDmThresholdMelting(double temperature_celsius, int species);
  static double getDmMeltingFactor(double dm, double dm_threshold,
                                   double transition_width);

  void meltingModelLiu24(double qr, double qs, double qg,
                         std::optional<double> ntr, std::optional<double> nts,
                         std::optional<double> ntg, std::optional<double> qh,
                         std::optional<double> nth,
                         const double qx_min[7], const double ntx_min[7],
                         std::optional<double> coeff_melt,
                         double density_air, double temperature,
                         MeltOut &o) const;

  void meltingModelJung08(double qr, double qs, double qg,
                          std::optional<double> ntr, std::optional<double> nts,
                          std::optional<double> ntg, std::optional<double> qh,
                          std::optional<double> nth,
                          MeltOut &o) const;
};

}  // namespace ppro
}  // namespace ufo

#endif  // UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_ZHANG21FORWARD_H_
