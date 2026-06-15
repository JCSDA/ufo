/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Author: Rong Kong (NCAR/MMM) - C++ rewrite of the operator from Fortran;
 *         bug fixes, operator tuning and testing, modularization and
 *         optimization.
 *
 * Original Fortran authors:
 *   Zhiquan (Jake) Liu (NCAR/MMM) - Initial implementation, July 2022.
 *   Hejun Xie (NCAR/MMM) - Integration into JEDI-UFO framework, August 2024.
 *   Tao Sun (NCAR/MMM) - Hail categories and extended microphysics, April 2025.
 */

#include "ufo/operators/radarpolarimetric/ppro/Zhang21Forward.h"

#include <algorithm>
#include <cmath>

#include "ufo/operators/radarpolarimetric/ppro/PPROCoefsData.h"

namespace ufo {
namespace ppro {

// Precip-type identifiers used by the WSM6 single-moment helpers.
enum { PT_RAIN = 0, PT_SNOW = 1, PT_GRAUPEL = 2, PT_HAIL = 3 };
// Species identifiers used by weticephaseDensity / getDmThresholdMelting.
enum { SP_SNOW = 0, SP_GRAUPEL = 1, SP_HAIL = 2 };

// ----------------------------------------------------------------------------
Zhang21Forward::Zhang21Forward(const PPROTuning &tuning) : t_(tuning) {
  // Precompute the dry-hydrometeor "a" coefficients.  With ratio = 0 (dry),
  // coef_a(c, 0) reduces to c(0) = iceCoef(coefs, i, 0).
  for (int i = 0; i < 13; ++i) {
    sband_snow_a_[i]    = iceCoef(sband_snow_coefs_data, i, 0);
    sband_graupel_a_[i] = iceCoef(sband_graupel_coefs_data, i, 0);
    sband_hail_a_[i]    = iceCoef(sband_hail_coefs_data, i, 0);
    cband_snow_a_[i]    = iceCoef(cband_snow_coefs_data, i, 0);
    cband_graupel_a_[i] = iceCoef(cband_graupel_coefs_data, i, 0);
    cband_hail_a_[i]    = iceCoef(cband_hail_coefs_data, i, 0);
  }
}

// ----------------------------------------------------------------------------
// Eq. (21) of Zhang et al., 2021
double Zhang21Forward::coefA(const double c[4], double gamma) {
  double a = c[0];
  double g = 1.0;
  for (int i = 1; i <= 3; ++i) {
    g *= gamma;
    a += c[i] * g;
  }
  return a;
}

// ----------------------------------------------------------------------------
// Eqs. (13) - (16) of Zhang et al., 2021
void Zhang21Forward::dualpolOpRain(double w, double dm, const double *a,
                                   double &zh, double &zdr, double &kdp,
                                   double &phv) {
  const double dm2 = dm * dm, dm3 = dm2 * dm, dm4 = dm3 * dm;
  const double poly1 = rainCoef(a, 1, 0) + rainCoef(a, 1, 1) * dm
                     + rainCoef(a, 1, 2) * dm2 + rainCoef(a, 1, 3) * dm3
                     + rainCoef(a, 1, 4) * dm4;
  zh = w * poly1 * poly1;

  zdr = rainCoef(a, 2, 0) + rainCoef(a, 2, 1) * dm + rainCoef(a, 2, 2) * dm2
      + rainCoef(a, 2, 3) * dm3 + rainCoef(a, 2, 4) * dm4;

  kdp = w * (rainCoef(a, 3, 0) + rainCoef(a, 3, 1) * dm + rainCoef(a, 3, 2) * dm2
           + rainCoef(a, 3, 3) * dm3 + rainCoef(a, 3, 4) * dm4);

  phv = rainCoef(a, 4, 0) + rainCoef(a, 4, 1) * dm + rainCoef(a, 4, 2) * dm2
      + rainCoef(a, 4, 3) * dm3 + rainCoef(a, 4, 4) * dm4;
}

// ----------------------------------------------------------------------------
// Eqs. (17) - (20) of Zhang et al., 2021 for snow/graupel/hail
void Zhang21Forward::dualpolOpIcephase(double z, double w, double rho, double dm,
                                       const double a[13],
                                       double &zh, double &zdr, double &kdp,
                                       double &phv, int iband,
                                       const std::string &ptype) {
  const double dm2 = dm * dm, dm3 = dm2 * dm;

  // C-band: the coefficient for C band cannot fit the S-band formula well for
  // hail due to strong non-Rayleigh scattering effects.
  if (iband == 2 && (ptype == "mh" || ptype == "ph")) {
    if (dm > 1.0e-3) {
      const double p = (a[0] + a[1] * dm + a[2] * dm2 + a[3] * dm3) / dm;
      zh = z * p * p;
    } else {
      zh = 0.0;
    }
  } else {
    const double p = a[0] + a[1] * dm + a[2] * dm2 + a[3] * dm3;
    zh = z * p * p;
  }
  zdr = a[4] + a[5] * dm + a[6] * dm2;

  if (rho > 1.0e-3) {
    kdp = w * (a[7] + a[8] * dm + a[9] * dm2) / rho;
  } else {
    kdp = 0.0;
  }
  phv = a[10] + a[11] * dm + a[12] * dm2;
}

// ----------------------------------------------------------------------------
// Eqs. 22 - 25 of Zhang et al., 2021
void Zhang21Forward::dualpolOpTotal(
    double zh_prain, double zh_psnow, double zh_pgraupel,
    double zdr_prain, double zdr_psnow, double zdr_pgraupel,
    double kdp_prain, double kdp_psnow, double kdp_pgraupel,
    double phv_prain, double phv_psnow, double phv_pgraupel,
    double &zh, double &zdr, double &kdp, double &phv,
    double zh_msnow, double zh_mgraupel,
    double zdr_msnow, double zdr_mgraupel,
    double kdp_msnow, double kdp_mgraupel,
    double phv_msnow, double phv_mgraupel,
    bool include_hail,
    double zh_mhail, double zdr_mhail, double kdp_mhail, double phv_mhail,
    double zh_phail, double zdr_phail, double kdp_phail, double phv_phail) {
  auto split = [](double zhx, double zdrx, double &zv, double &zbar) {
    if (zdrx >= 1.0) {
      zv = zhx / zdrx;
      zbar = zhx / std::sqrt(zdrx);
    } else {
      // Isotropic assumption (zdr < 1.0 non-physical here): zv = zh.
      zv = zhx;
      zbar = zhx;
    }
  };

  double zv_prain, zbar_prain, zv_msnow, zbar_msnow, zv_mgraupel, zbar_mgraupel;
  double zv_psnow, zbar_psnow, zv_pgraupel, zbar_pgraupel;
  split(zh_prain, zdr_prain, zv_prain, zbar_prain);
  split(zh_msnow, zdr_msnow, zv_msnow, zbar_msnow);
  split(zh_mgraupel, zdr_mgraupel, zv_mgraupel, zbar_mgraupel);
  split(zh_psnow, zdr_psnow, zv_psnow, zbar_psnow);
  split(zh_pgraupel, zdr_pgraupel, zv_pgraupel, zbar_pgraupel);

  double zv;
  if (include_hail) {
    double zv_mhail, zbar_mhail, zv_phail, zbar_phail;
    split(zh_mhail, zdr_mhail, zv_mhail, zbar_mhail);
    split(zh_phail, zdr_phail, zv_phail, zbar_phail);

    zh = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel
       + zh_mhail + zh_phail;
    zv = zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel
       + zv_mhail + zv_phail;
    zdr = (zv > 1.0e-3) ? zh / zv : 1.0;
    kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel
        + kdp_mhail + kdp_phail;
    const double zbar_sum = zbar_prain + zbar_msnow + zbar_mgraupel
                          + zbar_psnow + zbar_pgraupel + zbar_mhail + zbar_phail;
    if (zbar_sum > 1.0e-3) {
      phv = (zbar_prain * phv_prain + zbar_msnow * phv_msnow
           + zbar_mgraupel * phv_mgraupel + zbar_mhail * phv_mhail
           + zbar_psnow * phv_psnow + zbar_pgraupel * phv_pgraupel
           + zbar_phail * phv_phail) / zbar_sum;
    } else {
      phv = 1.0;
    }
  } else {
    zh = zh_prain + zh_msnow + zh_mgraupel + zh_psnow + zh_pgraupel;
    zv = zv_prain + zv_msnow + zv_mgraupel + zv_psnow + zv_pgraupel;
    zdr = (zv > 1.0e-3) ? zh / zv : 1.0;
    kdp = kdp_prain + kdp_msnow + kdp_mgraupel + kdp_psnow + kdp_pgraupel;
    const double zbar_sum = zbar_prain + zbar_msnow + zbar_mgraupel
                          + zbar_psnow + zbar_pgraupel;
    if (zbar_sum > 1.0e-3) {
      phv = (zbar_prain * phv_prain + zbar_msnow * phv_msnow
           + zbar_mgraupel * phv_mgraupel + zbar_psnow * phv_psnow
           + zbar_pgraupel * phv_pgraupel) / zbar_sum;
    } else {
      phv = 1.0;
    }
  }

  // KDP non-negative protection.
  kdp = std::max(kdp, 0.0);
}

// ----------------------------------------------------------------------------
// Eq. (7) of Zhang et al., 2021
double Zhang21Forward::weticephaseDensity(int species, double ratio) {
  const double r2 = ratio * ratio;
  switch (species) {
    case SP_SNOW:    return density_snow    * (1.0 - r2) + r2;
    case SP_GRAUPEL: return density_graupel * (1.0 - r2) + r2;
    case SP_HAIL:    return density_hail    * (1.0 - r2) + r2;
    default:         return density_snow    * (1.0 - r2) + r2;
  }
}

// ----------------------------------------------------------------------------
// Eqs. (5) and (6) of Zhang et al., 2021 (2-moment exponential PSD)
void Zhang21Forward::dmZ2moment(double air_density, double qx, double rhox,
                                double ntx, double &dm, double &z,
                                std::optional<double> dmmax,
                                std::optional<double> dm_min) const {
  const double thrd = 1.0 / 3.0;
  if (t_.enable_dm_regularization) {
    if (qx < 1.0e-3 && ntx < 1000.0) {
      dm = 0.0;
    } else {
      dm = 4.0 * std::pow((air_density * qx) / (M_PI * rhox * ntx), thrd);
      dm = dm * 10.0;  // cm -> mm
    }
  } else {
    if (ntx < 1.0) {
      dm = 0.0;
    } else {
      dm = 4.0 * std::pow((air_density * qx) / (M_PI * rhox * ntx), thrd);
      dm = dm * 10.0;  // cm -> mm
    }
  }

  if (dmmax) dm = std::min(dm, *dmmax);
  if (dm_min && dm < *dm_min) dm = 0.0;

  z = 11250.0 * ((air_density * qx) / (M_PI * rhox)) * dm * dm * dm;
}

// ----------------------------------------------------------------------------
// Single-moment exponential PSD: N(D) = N0 exp(-lambda D)
void Zhang21Forward::n0LambdaWsm6(int precip_type, double air_density,
                                  double temp, double qx, double rhox,
                                  double &n0, double &lambda) {
  switch (precip_type) {
    case PT_RAIN:
      n0 = 8.0e6;
      break;
    case PT_SNOW:
      n0 = 2.0e6 * std::exp(0.12 * (273.15 - temp));  // Hong et al. 2004 Eq.(6)
      break;
    case PT_GRAUPEL:
    case PT_HAIL:
      n0 = 4.0e6;
      break;
    default:
      n0 = 4.0e6;
      break;
  }
  lambda = std::pow((M_PI * 1000.0 * rhox * n0) / (air_density * qx * 0.001), 0.25);
}

// ----------------------------------------------------------------------------
void Zhang21Forward::dmZwsm6(int precip_type, double air_density, double temp,
                             double qx, double rhox, double &dm, double &z,
                             std::optional<double> dmmax,
                             std::optional<double> dm_min) {
  // Guard against qx <= 0: n0LambdaWsm6 divides by qx to obtain lambda.
  // A non-positive mass yields no signal, matching the qx < 1.0e-3 cutoff below.
  if (qx <= 0.0) {
    dm = 0.0;
    z = 0.0;
    return;
  }
  double n0, lambda;
  n0LambdaWsm6(precip_type, air_density, temp, qx, rhox, n0, lambda);

  dm = 4000.0 / lambda;  // Eq. (5), mm

  if (dmmax) dm = std::min(dm, *dmmax);
  if (dm_min && dm < *dm_min) dm = 0.0;

  z = 11250.0 * ((air_density * qx) / (M_PI * rhox)) * dm * dm * dm;  // Eq. (6)

  if (qx < 1.0e-3) {
    dm = 0.0;
    z = 0.0;
  }
}

// ----------------------------------------------------------------------------
double Zhang21Forward::getDmThresholdMelting(double tc, int species) {
  const double thresh1 = -5.0;  // cold-end T [C]
  const double thresh2 =  0.0;  // warm-end T [C]
  double dm_cold, dm_warm;
  switch (species) {
    case SP_SNOW:    dm_cold = 1.0; dm_warm = 2.0; break;
    case SP_GRAUPEL: dm_cold = 1.5; dm_warm = 3.0; break;
    case SP_HAIL:    dm_cold = 2.0; dm_warm = 4.0; break;
    default:         dm_cold = 1.0; dm_warm = 2.0; break;
  }
  if (tc < thresh1) return dm_cold;
  if (tc > thresh2) return dm_warm;
  return dm_cold + (tc - thresh1) / (thresh2 - thresh1) * (dm_warm - dm_cold);
}

// ----------------------------------------------------------------------------
double Zhang21Forward::getDmMeltingFactor(double dm, double dm_threshold,
                                          double transition_width) {
  if (transition_width > 0.0) {
    return 1.0 / (1.0 + std::exp((dm - dm_threshold) / transition_width));
  }
  return (dm < dm_threshold) ? 1.0 : 0.0;
}

// ----------------------------------------------------------------------------
// Melting model based on Liu et al., 2024 with modifications by Kong (NCAR/MMM)
void Zhang21Forward::meltingModelLiu24(
    double qr, double qs, double qg,
    std::optional<double> ntr, std::optional<double> nts,
    std::optional<double> ntg, std::optional<double> qh,
    std::optional<double> nth,
    const double qx_min[7], const double ntx_min[7],
    std::optional<double> coeff_melt,
    double density_air, double temperature, MeltOut &o) const {
  const bool double_moment = ntr || nts || ntg || nth;

  double factor_s = 1.0, factor_g = 1.0, factor_h = 1.0;
  double eff_coeff_s = 0.0, eff_coeff_g = 0.0, eff_coeff_h = 0.0;
  const double melt_coeff = coeff_melt.value_or(0.3);

  double dm_factor_s = 1.0, dm_factor_g = 1.0, dm_factor_h = 1.0;

  if (t_.enable_dm_melting_limit) {
    const double dm_thresh_snow = getDmThresholdMelting(temperature, SP_SNOW);
    const double dm_thresh_graupel = getDmThresholdMelting(temperature, SP_GRAUPEL);
    const double dm_thresh_hail = getDmThresholdMelting(temperature, SP_HAIL);
    double dmp, ztmp;

    if (qs > 1.0e-6) {
      if (nts && *nts > 1.0e-6) {
        dmZ2moment(density_air, qs, density_snow, *nts, dmp, ztmp,
                   t_.dmmax_pure_snow, std::nullopt);
      } else {
        dmZwsm6(PT_SNOW, density_air, temperature + 273.15, qs, density_snow,
                dmp, ztmp, t_.dmmax_pure_snow, std::nullopt);
      }
      dm_factor_s = (dmp >= 0.3)
          ? getDmMeltingFactor(dmp, dm_thresh_snow, t_.dm_melting_transition_width)
          : 0.0;
    }
    if (qg > 1.0e-6) {
      if (ntg && *ntg > 1.0e-6) {
        dmZ2moment(density_air, qg, density_graupel, *ntg, dmp, ztmp,
                   t_.dmmax_pure_graupel, std::nullopt);
      } else {
        dmZwsm6(PT_GRAUPEL, density_air, temperature + 273.15, qg, density_graupel,
                dmp, ztmp, t_.dmmax_pure_graupel, std::nullopt);
      }
      dm_factor_g = (dmp >= 0.3)
          ? getDmMeltingFactor(dmp, dm_thresh_graupel, t_.dm_melting_transition_width)
          : 0.0;
    }
    if (qh && *qh > 1.0e-6) {
      if (nth && *nth > 1.0e-6) {
        dmZ2moment(density_air, *qh, density_hail, *nth, dmp, ztmp,
                   t_.dmmax_pure_hail, std::nullopt);
      } else {
        dmZwsm6(PT_HAIL, density_air, temperature + 273.15, *qh, density_hail,
                dmp, ztmp, t_.dmmax_pure_hail, std::nullopt);
      }
      dm_factor_h = (dmp >= 0.3)
          ? getDmMeltingFactor(dmp, dm_thresh_hail, t_.dm_melting_transition_width)
          : 0.0;
    }
  }

  // Snow melting
  if (qr + qs < qsum_min_threshold) {
    o.qms = 0.0;
    o.rats = 0.0;
  } else {
    if (t_.enable_melting_transition) {
      // Guard the ratio: if either category is zero the species are maximally
      // unbalanced (balance = 0). Avoids division by zero under FPE trapping.
      const double balance = (qr > 0.0 && qs > 0.0)
                           ? std::min(qs / qr, qr / qs) : 0.0;
      if (balance < t_.snow_ratio_low) {
        factor_s = 0.0;
      } else if (balance > t_.snow_ratio_high) {
        factor_s = 1.0;
      } else if (t_.snow_ratio_high <= t_.snow_ratio_low) {
        factor_s = 1.0;
      } else {
        factor_s = (balance - t_.snow_ratio_low)
                 / (t_.snow_ratio_high - t_.snow_ratio_low);
      }
    }
    eff_coeff_s = melt_coeff * factor_s * dm_factor_s;
    o.qms = std::pow(qr, t_.melting_exp_rain) * std::pow(qs, t_.melting_exp_ice)
          * eff_coeff_s;
    o.rats = std::pow(qr, t_.fraction_exp_rain)
           / (std::pow(qr, t_.fraction_exp_rain) + std::pow(qs, t_.fraction_exp_ice));
    o.qms = std::max(o.qms, qx_min[1]);
  }

  // Graupel melting
  if (qr + qg < qsum_min_threshold) {
    o.qmg = 0.0;
    o.ratg = 0.0;
  } else {
    if (t_.enable_melting_transition) {
      // Guard the ratio: if either category is zero the species are maximally
      // unbalanced (balance = 0). Avoids division by zero under FPE trapping.
      const double balance = (qr > 0.0 && qg > 0.0)
                           ? std::min(qg / qr, qr / qg) : 0.0;
      if (balance < t_.graupel_ratio_low) {
        factor_g = 0.0;
      } else if (balance > t_.graupel_ratio_high) {
        factor_g = 1.0;
      } else if (t_.graupel_ratio_high <= t_.graupel_ratio_low) {
        factor_g = 1.0;
      } else {
        factor_g = (balance - t_.graupel_ratio_low)
                 / (t_.graupel_ratio_high - t_.graupel_ratio_low);
      }
    }
    eff_coeff_g = melt_coeff * factor_g * dm_factor_g;
    o.qmg = std::pow(qr, t_.melting_exp_rain) * std::pow(qg, t_.melting_exp_ice)
          * eff_coeff_g;
    o.ratg = std::pow(qr, t_.fraction_exp_rain)
           / (std::pow(qr, t_.fraction_exp_rain) + std::pow(qg, t_.fraction_exp_ice));
    o.qmg = std::max(o.qmg, qx_min[2]);
  }

  // Hail melting
  if (qh) {
    if (qr + *qh < qsum_min_threshold) {
      o.qmh = 0.0;
      o.rath = 0.0;
    } else {
      if (t_.enable_melting_transition) {
        // Guard the ratio: if either category is zero the species are maximally
        // unbalanced (balance = 0). Avoids division by zero under FPE trapping.
        const double balance = (qr > 0.0 && *qh > 0.0)
                             ? std::min(*qh / qr, qr / *qh) : 0.0;
        if (balance < t_.hail_ratio_low) {
          factor_h = 0.0;
        } else if (balance > t_.hail_ratio_high) {
          factor_h = 1.0;
        } else if (t_.hail_ratio_high <= t_.hail_ratio_low) {
          factor_h = 1.0;
        } else {
          factor_h = (balance - t_.hail_ratio_low)
                   / (t_.hail_ratio_high - t_.hail_ratio_low);
        }
      }
      eff_coeff_h = melt_coeff * factor_h * dm_factor_h;
      o.qmh = std::pow(qr, t_.melting_exp_rain) * std::pow(*qh, t_.melting_exp_ice)
            * eff_coeff_h;
      o.rath = std::pow(qr, t_.fraction_exp_rain)
             / (std::pow(qr, t_.fraction_exp_rain) + std::pow(*qh, t_.fraction_exp_ice));
      o.qmh = std::max(o.qmh, qx_min[3]);
    }
  }

  // Scale melting fractions by melting_water_fraction (default 1.0).
  o.rats = std::min(t_.melting_water_fraction * o.rats, 1.0);
  o.ratg = std::min(t_.melting_water_fraction * o.ratg, 1.0);
  if (qh) o.rath = std::min(t_.melting_water_fraction * o.rath, 1.0);

  const double qmsr = o.qms * o.rats;
  const double qmgr = o.qmg * o.ratg;
  double qmhr = 0.0;
  if (qh) qmhr = o.qmh * o.rath;

  const double qmsd = o.qms - qmsr;
  const double qmgd = o.qmg - qmgr;
  double qmhd = 0.0;
  if (qh) qmhd = o.qmh - qmhr;

  o.qpr = qr - qmsr - qmgr;
  o.qps = qs - qmsd;
  o.qpg = qg - qmgd;
  if (qh) {
    o.qpr = o.qpr - qmhr;
    o.qph = *qh - qmhd;
  }

  o.qpr = std::max(o.qpr, qx_min[0]);
  o.qps = std::max(o.qps, qx_min[4]);
  o.qpg = std::max(o.qpg, qx_min[5]);
  if (qh) o.qph = std::max(o.qph, qx_min[6]);

  if (double_moment) {
    if (nts && ntr) {
      if (o.qms > 1.0e-3) {
        o.ntms = std::pow(*ntr, t_.melting_exp_rain) * std::pow(*nts, t_.melting_exp_ice)
               * eff_coeff_s;
        o.ntms = std::max(o.ntms, ntx_min[1]);
      } else {
        o.ntms = ntx_min[1];
      }
    }
    if (ntg && ntr) {
      if (o.qmg > 1.0e-3) {
        o.ntmg = std::pow(*ntr, t_.melting_exp_rain) * std::pow(*ntg, t_.melting_exp_ice)
               * eff_coeff_g;
        o.ntmg = std::max(o.ntmg, ntx_min[2]);
      } else {
        o.ntmg = ntx_min[2];
      }
    }
    if (nth && ntr) {
      if (qh && o.qmh > 1.0e-3) {
        o.ntmh = std::pow(*ntr, t_.melting_exp_rain) * std::pow(*nth, t_.melting_exp_ice)
               * eff_coeff_h;
        o.ntmh = std::max(o.ntmh, ntx_min[3]);
      } else {
        o.ntmh = ntx_min[3];
      }
    }
    if (ntr) o.ntpr = std::max(*ntr, ntx_min[0]);
    if (nts) o.ntps = std::max(*nts, ntx_min[4]);
    if (ntg) o.ntpg = std::max(*ntg, ntx_min[5]);
    if (nth) o.ntph = std::max(*nth, ntx_min[6]);
  }
}

// ----------------------------------------------------------------------------
// Classic melting scheme based on Jung et al., 2008
void Zhang21Forward::meltingModelJung08(
    double qr, double qs, double qg,
    std::optional<double> ntr, std::optional<double> nts,
    std::optional<double> ntg, std::optional<double> qh,
    std::optional<double> nth, MeltOut &o) const {
  const bool double_moment = ntr && nts && ntg;

  if (qs > 0.0 && qr > 0.0) {
    const double minratio = std::min(qr / qs, qs / qr);
    if (minratio > 0.01) {
      o.qms = (qr + qs) * 0.5 * std::pow(minratio, 0.3);
      o.rats = qr / (qr + qs);
    } else {
      o.qms = 0.0;
      o.rats = 0.0;
    }
  } else {
    o.qms = 0.0;
    o.rats = 0.0;
  }

  if (qg > 0.0 && qr > 0.0) {
    const double minratio = std::min(qr / qg, qg / qr);
    if (minratio > 0.01) {
      o.qmg = (qr + qg) * 0.5 * std::pow(minratio, 0.3);
      o.ratg = qr / (qr + qg);
    } else {
      o.qmg = 0.0;
      o.ratg = 0.0;
    }
  } else {
    o.qmg = 0.0;
    o.ratg = 0.0;
  }

  if (qh) {
    if (*qh > 0.0 && qr > 0.0) {
      const double minratio = std::min(qr / *qh, *qh / qr);
      if (minratio > 0.01) {
        o.qmh = (qr + *qh) * 0.5 * std::pow(minratio, 0.3);
        o.rath = qr / (qr + *qh);
      } else {
        o.qmh = 0.0;
        o.rath = 0.0;
      }
    }
  }

  const double qmsr = o.qms * o.rats;
  const double qmgr = o.qmg * o.ratg;
  const double qmsd = o.qms * (1.0 - o.rats);
  const double qmgd = o.qmg * (1.0 - o.ratg);
  double qmhr = 0.0, qmhd = 0.0;
  if (qh) {
    qmhr = o.qmh * o.rath;
    qmhd = o.qmh * (1.0 - o.rath);
  }

  o.qpr = qr - qmsr - qmgr;
  o.qps = qs - qmsd;
  o.qpg = qg - qmgd;
  if (qh) {
    o.qpr = o.qpr - qmhr;
    o.qph = *qh - qmhd;
  }

  if (double_moment) {
    o.ntms = std::sqrt(*ntr * *nts);
    o.ntmg = std::sqrt(*ntr * *ntg);
    o.ntpr = *ntr;
    o.ntps = *nts;
    o.ntpg = *ntg;
    if (nth) {
      o.ntmh = std::sqrt(*ntr * *nth);
      o.ntph = *nth;
    }
  }
}

// ============================================================================
// High-level interface: compute dual-pol variables for one observation point.
void Zhang21Forward::computePoint(int iband, double density_air, double temp_air,
                                  double qr, double qs, double qg,
                                  double &zh, double &zdr, double &kdp,
                                  double &phv,
                                  const Zhang21Optionals &opt) const {
  const double temperature_celsius =
      opt.temperature ? *opt.temperature : (temp_air - 273.15);

  double qx_min[7];
  double ntx_min[7];
  for (int i = 0; i < 7; ++i) {
    qx_min[i] = 0.0;
    ntx_min[i] = 10.0;
  }
  // Index: 1=pr, 2=ms, 3=mg, 4=mh, 5=ps, 6=pg, 7=ph (1-based -> 0-based here).
  double dmmax[7];
  dmmax[0] = t_.dmmax_rain;
  dmmax[1] = t_.dmmax_melting_snow;
  dmmax[2] = t_.dmmax_melting_graupel;
  dmmax[3] = t_.dmmax_melting_hail;
  dmmax[4] = t_.dmmax_pure_snow;
  dmmax[5] = t_.dmmax_pure_graupel;
  dmmax[6] = t_.dmmax_pure_hail;

  const double *rain_coefs, *snow_coefs, *graupel_coefs, *hail_coefs;
  const double *snow_a, *graupel_a, *hail_a;
  if (iband == 1) {
    rain_coefs = sband_rain_coefs_data;
    snow_coefs = sband_snow_coefs_data;
    graupel_coefs = sband_graupel_coefs_data;
    hail_coefs = sband_hail_coefs_data;
    snow_a = sband_snow_a_;
    graupel_a = sband_graupel_a_;
    hail_a = sband_hail_a_;
  } else {
    rain_coefs = cband_rain_coefs_data;
    snow_coefs = cband_snow_coefs_data;
    graupel_coefs = cband_graupel_coefs_data;
    hail_coefs = cband_hail_coefs_data;
    snow_a = cband_snow_a_;
    graupel_a = cband_graupel_a_;
    hail_a = cband_hail_a_;
  }

  // Regularize inputs.
  const double qrreg = std::max(qx_min[0], qr);
  const double qsreg = std::max(qx_min[4], qs);
  const double qgreg = std::max(qx_min[5], qg);
  std::optional<double> qhreg, nrreg, nsreg, ngreg, nhreg;
  if (opt.qh) qhreg = std::max(qx_min[6], *opt.qh);
  if (opt.nr) nrreg = std::max(ntx_min[0], *opt.nr);
  if (opt.ns) nsreg = std::max(ntx_min[4], *opt.ns);
  if (opt.ng) ngreg = std::max(ntx_min[5], *opt.ng);
  if (opt.nh) nhreg = std::max(ntx_min[6], *opt.nh);

  // Dispatch melting model.
  MeltOut m;
  if (t_.melting_scheme_option == "jung08") {
    meltingModelJung08(qrreg, qsreg, qgreg, nrreg, nsreg, ngreg, qhreg, nhreg, m);
  } else {
    meltingModelLiu24(qrreg, qsreg, qgreg, nrreg, nsreg, ngreg, qhreg, nhreg,
                      qx_min, ntx_min, opt.coeff_melt, density_air,
                      temperature_celsius, m);
  }

  // Apply melting-fraction tuning coefficients (downstream scattering only).
  if (t_.tuning_melt_frac_snow != 1.0)
    m.rats = std::min(t_.tuning_melt_frac_snow * m.rats, 1.0);
  if (t_.tuning_melt_frac_graupel != 1.0)
    m.ratg = std::min(t_.tuning_melt_frac_graupel * m.ratg, 1.0);
  if (opt.qh && t_.tuning_melt_frac_hail != 1.0)
    m.rath = std::min(t_.tuning_melt_frac_hail * m.rath, 1.0);

  // Regularize melting-scheme outputs.
  m.qpr = std::max(qx_min[0], m.qpr);
  m.qps = std::max(qx_min[4], m.qps);
  m.qpg = std::max(qx_min[5], m.qpg);
  if (opt.qh && opt.nh) {
    m.qph = std::max(qx_min[6], m.qph);
    m.ntpr = std::max(ntx_min[0], m.ntpr);
    m.ntms = std::max(ntx_min[1], m.ntms);
    m.ntmg = std::max(ntx_min[2], m.ntmg);
    m.ntmh = std::max(ntx_min[3], m.ntmh);
    m.ntps = std::max(ntx_min[4], m.ntps);
    m.ntpg = std::max(ntx_min[5], m.ntpg);
    m.ntph = std::max(ntx_min[6], m.ntph);
  } else if (opt.nr && opt.ns && opt.ng) {
    m.ntpr = std::max(ntx_min[0], m.ntpr);
    m.ntms = std::max(ntx_min[1], m.ntms);
    m.ntmg = std::max(ntx_min[2], m.ntmg);
    m.ntps = std::max(ntx_min[4], m.ntps);
    m.ntpg = std::max(ntx_min[5], m.ntpg);
  } else if (opt.nr) {
    m.ntpr = std::max(ntx_min[0], m.ntpr);
  }

  double zhpr, zdrpr, kdppr, phvpr;
  double zhms, zdrms, kdpms, phvms;
  double zhps, zdrps, kdpps, phvps;
  double zhmg, zdrmg, kdpmg, phvmg;
  double zhpg, zdrpg, kdppg, phvpg;
  double zhmh = 0.0, zdrmh = 1.0, kdpmh = 0.0, phvmh = 1.0;
  double zhph = 0.0, zdrph = 1.0, kdpph = 0.0, phvph = 1.0;
  double wpr, zpr, dmpr, wms, zms, dmms, denms;
  double wps, zps, dmps, wmg, zmg, dmmg, denmg, wpg, zpg, dmpg, denpg;
  double c4[4];

  // 1. Pure rain
  if (t_.skip_small_qx && m.qpr < 1.0e-3) {
    wpr = 0.0; zpr = 0.0; dmpr = 0.0;
    zhpr = 0.0; zdrpr = 1.0; kdppr = 0.0; phvpr = 1.0;
  } else {
    wpr = waterContent(density_air, m.qpr);
    if (opt.nr) {
      dmZ2moment(density_air, m.qpr, density_rain, m.ntpr, dmpr, zpr,
                 dmmax[0], std::nullopt);
    } else {
      dmZwsm6(PT_RAIN, density_air, temp_air, m.qpr, density_rain, dmpr, zpr,
              dmmax[0], std::nullopt);
    }
    if (t_.tuning_dm_rain != 1.0) {
      dmpr = t_.tuning_dm_rain * dmpr;
      zpr = zpr * t_.tuning_dm_rain * t_.tuning_dm_rain * t_.tuning_dm_rain;
    }
    dualpolOpRain(wpr, dmpr, rain_coefs, zhpr, zdrpr, kdppr, phvpr);
  }

  // 2. Melting snow
  if (t_.skip_small_qx && m.qms < 0.001) {
    wms = 0.0; zms = 0.0; dmms = 0.0;
    denms = density_snow;
    zhms = 0.0; zdrms = 1.0; kdpms = 0.0; phvms = 1.0;
  } else {
    wms = waterContent(density_air, m.qms);
    denms = weticephaseDensity(SP_SNOW, m.rats);
    denms = std::min(0.5, std::max(0.1, denms));
    if (opt.ns) {
      dmZ2moment(density_air, m.qms, denms, m.ntms, dmms, zms, dmmax[1], 0.3);
    } else {
      dmZwsm6(PT_SNOW, density_air, temp_air, m.qms, denms, dmms, zms, dmmax[1], 0.3);
    }
    if (t_.tuning_dm_melting_snow != 1.0) {
      dmms = t_.tuning_dm_melting_snow * dmms;
      zms = zms * t_.tuning_dm_melting_snow * t_.tuning_dm_melting_snow
          * t_.tuning_dm_melting_snow;
    }
    double as[13];
    for (int i = 0; i < 13; ++i) {
      for (int j = 0; j < 4; ++j) c4[j] = iceCoef(snow_coefs, i, j);
      as[i] = coefA(c4, m.rats);
    }
    dualpolOpIcephase(zms, wms, denms, dmms, as, zhms, zdrms, kdpms, phvms, 0, "");
  }

  // 3. Pure snow
  if (t_.skip_small_qx && m.qps < 1.0e-4) {
    wps = 0.0; zps = 0.0; dmps = 0.0;
    zhps = 0.0; zdrps = 1.0; kdpps = 0.0; phvps = 1.0;
  } else {
    wps = waterContent(density_air, m.qps);
    if (opt.ns) {
      dmZ2moment(density_air, m.qps, density_snow, m.ntps, dmps, zps, dmmax[4], 0.3);
    } else {
      dmZwsm6(PT_SNOW, density_air, temp_air, m.qps, density_snow, dmps, zps,
              dmmax[4], 0.3);
    }
    if (t_.tuning_dm_pure_snow != 1.0) {
      dmps = t_.tuning_dm_pure_snow * dmps;
      zps = zps * t_.tuning_dm_pure_snow * t_.tuning_dm_pure_snow
          * t_.tuning_dm_pure_snow;
    }
    dualpolOpIcephase(zps, wps, density_snow, dmps, snow_a, zhps, zdrps, kdpps,
                      phvps, 0, "");
  }

  // 4. Melting graupel
  if (t_.skip_small_qx && m.qmg < 0.01) {
    wmg = 0.0; zmg = 0.0; dmmg = 0.0;
    denmg = density_graupel;
    zhmg = 0.0; zdrmg = 1.0; kdpmg = 0.0; phvmg = 1.0;
  } else {
    wmg = waterContent(density_air, m.qmg);
    denmg = weticephaseDensity(SP_GRAUPEL, m.ratg);
    if (opt.ng) {
      dmZ2moment(density_air, m.qmg, denmg, m.ntmg, dmmg, zmg, dmmax[2], 0.3);
    } else {
      dmZwsm6(PT_GRAUPEL, density_air, temp_air, m.qmg, denmg, dmmg, zmg, dmmax[2], 0.3);
    }
    if (t_.tuning_dm_melting_graupel != 1.0) {
      dmmg = t_.tuning_dm_melting_graupel * dmmg;
      zmg = zmg * t_.tuning_dm_melting_graupel * t_.tuning_dm_melting_graupel
          * t_.tuning_dm_melting_graupel;
    }
    double ag[13];
    for (int i = 0; i < 13; ++i) {
      for (int j = 0; j < 4; ++j) c4[j] = iceCoef(graupel_coefs, i, j);
      ag[i] = coefA(c4, m.ratg);
    }
    dualpolOpIcephase(zmg, wmg, denmg, dmmg, ag, zhmg, zdrmg, kdpmg, phvmg, 0, "");
  }

  // 5. Pure graupel
  if (t_.skip_small_qx && m.qpg < 0.001) {
    wpg = 0.0; zpg = 0.0; dmpg = 0.0;
    denpg = density_graupel;
    zhpg = 0.0; zdrpg = 1.0; kdppg = 0.0; phvpg = 1.0;
  } else {
    wpg = waterContent(density_air, m.qpg);
    if (opt.vg && opt.ng) {
      if (m.qpg > 0.0 && *opt.vg > 0.0 && *opt.ng > 0.0) {
        denpg = density_air * m.qpg / *opt.vg * 1.0e-6;
        denpg = std::min(std::max(denpg, 0.2), 0.6);
      } else {
        denpg = density_graupel;
      }
    } else {
      denpg = density_graupel;
    }
    if (opt.ng) {
      dmZ2moment(density_air, m.qpg, denpg, m.ntpg, dmpg, zpg, dmmax[5], 0.3);
    } else {
      dmZwsm6(PT_GRAUPEL, density_air, temp_air, m.qpg, denpg, dmpg, zpg, dmmax[5], 0.3);
    }
    if (t_.tuning_dm_pure_graupel != 1.0) {
      dmpg = t_.tuning_dm_pure_graupel * dmpg;
      zpg = zpg * t_.tuning_dm_pure_graupel * t_.tuning_dm_pure_graupel
          * t_.tuning_dm_pure_graupel;
    }
    dualpolOpIcephase(zpg, wpg, denpg, dmpg, graupel_a, zhpg, zdrpg, kdppg,
                      phvpg, 0, "");
  }

  const bool include_hail = opt.nh && opt.qh;

  if (include_hail) {
    double wmh, zmh, dmmh, denmh, wph, zph, dmph, denph;

    // 6. Melting hail
    if (t_.skip_small_qx && m.qmh < 0.01) {
      wmh = 0.0; zmh = 0.0; dmmh = 0.0;
      denmh = t_.treat_hail_as_graupel ? density_graupel : density_hail;
      zhmh = 0.0; zdrmh = 1.0; kdpmh = 0.0; phvmh = 1.0;
    } else {
      wmh = waterContent(density_air, m.qmh);
      if (t_.treat_hail_as_graupel) {
        denmh = weticephaseDensity(SP_GRAUPEL, m.rath);
      } else {
        denmh = weticephaseDensity(SP_HAIL, m.rath);
      }
      if (opt.nh) {
        dmZ2moment(density_air, m.qmh, denmh, m.ntmh, dmmh, zmh, dmmax[3], 0.3);
      } else if (t_.treat_hail_as_graupel) {
        dmZwsm6(PT_GRAUPEL, density_air, temp_air, m.qmh, denmh, dmmh, zmh, dmmax[3], 0.3);
      } else {
        dmZwsm6(PT_HAIL, density_air, temp_air, m.qmh, denmh, dmmh, zmh, dmmax[3], 0.3);
      }
      if (t_.tuning_dm_melting_hail != 1.0) {
        dmmh = t_.tuning_dm_melting_hail * dmmh;
        zmh = zmh * t_.tuning_dm_melting_hail * t_.tuning_dm_melting_hail
            * t_.tuning_dm_melting_hail;
      }
      double ah[13];
      for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 4; ++j) {
          c4[j] = t_.treat_hail_as_graupel ? iceCoef(graupel_coefs, i, j)
                                           : iceCoef(hail_coefs, i, j);
        }
        ah[i] = coefA(c4, m.rath);
      }
      if (t_.treat_hail_as_graupel) {
        dualpolOpIcephase(zmh, wmh, denmh, dmmh, ah, zhmh, zdrmh, kdpmh, phvmh,
                          iband, "");
      } else {
        dualpolOpIcephase(zmh, wmh, denmh, dmmh, ah, zhmh, zdrmh, kdpmh, phvmh,
                          iband, "mh");
      }
    }

    // 7. Pure hail
    if (t_.skip_small_qx && m.qph < 0.001) {
      wph = 0.0; zph = 0.0; dmph = 0.0;
      denph = t_.treat_hail_as_graupel ? density_graupel : density_hail;
      zhph = 0.0; zdrph = 1.0; kdpph = 0.0; phvph = 1.0;
    } else {
      wph = waterContent(density_air, m.qph);
      if (opt.vh && opt.nh) {
        if (m.qph > 0.0 && *opt.vh > 0.0 && *opt.nh > 0.0) {
          denph = density_air * m.qph / *opt.vh * 1.0e-6;
          denph = std::min(std::max(denph, 0.7), 1.0);
        } else {
          denph = t_.treat_hail_as_graupel ? density_graupel : density_hail;
        }
      } else {
        denph = t_.treat_hail_as_graupel ? density_graupel : density_hail;
      }
      if (opt.nh) {
        dmZ2moment(density_air, m.qph, denph, m.ntph, dmph, zph, dmmax[6], 0.3);
      } else if (t_.treat_hail_as_graupel) {
        dmZwsm6(PT_GRAUPEL, density_air, temp_air, m.qph, denph, dmph, zph, dmmax[6], 0.3);
      } else {
        dmZwsm6(PT_HAIL, density_air, temp_air, m.qph, denph, dmph, zph, dmmax[6], 0.3);
      }
      if (t_.tuning_dm_pure_hail != 1.0) {
        dmph = t_.tuning_dm_pure_hail * dmph;
        zph = zph * t_.tuning_dm_pure_hail * t_.tuning_dm_pure_hail
            * t_.tuning_dm_pure_hail;
      }
      if (t_.treat_hail_as_graupel) {
        double ah[13];
        for (int i = 0; i < 13; ++i) {
          for (int j = 0; j < 4; ++j) c4[j] = iceCoef(graupel_coefs, i, j);
          ah[i] = coefA(c4, 0.0);
        }
        dualpolOpIcephase(zph, wph, denph, dmph, ah, zhph, zdrph, kdpph, phvph,
                          iband, "");
      } else {
        dualpolOpIcephase(zph, wph, denph, dmph, hail_a, zhph, zdrph, kdpph,
                          phvph, iband, "ph");
      }
    }

    dualpolOpTotal(zhpr, zhps, zhpg, zdrpr, zdrps, zdrpg,
                   kdppr, kdpps, kdppg, phvpr, phvps, phvpg,
                   zh, zdr, kdp, phv,
                   zhms, zhmg, zdrms, zdrmg, kdpms, kdpmg, phvms, phvmg,
                   true,
                   zhmh, zdrmh, kdpmh, phvmh, zhph, zdrph, kdpph, phvph);
  } else {
    dualpolOpTotal(zhpr, zhps, zhpg, zdrpr, zdrps, zdrpg,
                   kdppr, kdpps, kdppg, phvpr, phvps, phvpg,
                   zh, zdr, kdp, phv,
                   zhms, zhmg, zdrms, zdrmg, kdpms, kdpmg, phvms, phvmg,
                   false,
                   0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  }
}

}  // namespace ppro
}  // namespace ufo
